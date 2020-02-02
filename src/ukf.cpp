#include "ukf.h"
#include "Eigen/Dense"
#include "iostream"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  if(this->is_initialized_)
  {
    return;
  }
  std::cout <<"Enter UKF constructor\n";
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  this->Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  this->weights_ = VectorXd(2 * n_aug_ + 1);
  this->is_initialized_ = false;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

	/*****************************************************************************
	*  Initialization
	****************************************************************************/
	if (!is_initialized_) {
		// first measurement
		this->x_ << 1, 1, 1, 1, 1;//Random init values for state vector.

		float x_reading = 1, y_reading = 1;
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
			const float ro = meas_package.raw_measurements_(0);
			const float theta = meas_package.raw_measurements_(1);
			float ro_dot = meas_package.raw_measurements_(2);
			x_reading = ro * cos(theta);
			y_reading = ro * sin(theta);
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			/**
			Initialize state.
			*/
			x_reading = meas_package.raw_measurements_(0);
			y_reading = meas_package.raw_measurements_(1);
		}
		this->x_(0) = x_reading;
		this->x_(1) = y_reading;

		P_ << 1, 0, 0, 0, 0,
			0, 1, 0, 0, 0,
			0, 0, 1, 0, 0,
			0, 0, 0, 1, 0,
			0, 0, 0, 0, 1;
		this->previous_timestamp_ = meas_package.timestamp_;

		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}
	/*****************************************************************************
	*  Prediction
	****************************************************************************/

	/**
	TODO:
	* Update the state transition matrix F according to the new elapsed time.
	- Time is measured in seconds.
	* Update the process noise covariance matrix.
	* Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
	*/
	//compute the time elapsed between the current and previous measurements
	float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = meas_package.timestamp_;

	this->Prediction(dt);

	/*****************************************************************************
	*  Update
	****************************************************************************/

	/**
	TODO:
	* Use the sensor type to perform the update step.
	* Update the state and covariance matrices.
	*/

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		// Radar updates
		this->UpdateRadar(meas_package);
	}
	else {
		// Laser updates
		this->UpdateLidar(meas_package);
	}

	// print the output
	//std::cout << "x_ = " << this->x_ << endl;
	//std::cout << "P_ = " << this->P_ << endl;
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  /*******************************************************************************
  * Generate segma points
  ******************************************************************************/
  // Lesson 7. section 18
  //define spreading parameter
	double lambda = 3 - this->n_x_;

	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	//create augmented mean vector
	VectorXd x_aug = VectorXd(7);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(7, 7);

	//create augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5, 5) = P_;
	P_aug(5, 5) = std_a_ * std_a_;
	P_aug(6, 6) = std_yawdd_ * std_yawdd_;

	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i< n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda + n_aug_) * L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda + n_aug_) * L.col(i);
	}

	//print result
	//std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

	/*******************************************************************************
	* Predict segma points
	******************************************************************************/

	//Lesson 7, section 21.

	//predict sigma points
	for (int i = 0; i< 2 * n_aug_ + 1; i++)
	{
		//extract values for better readability
		double p_x = Xsig_aug(0, i);
		double p_y = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
		}
		else {
			px_p = p_x + v * delta_t * cos(yaw);
			py_p = p_y + v * delta_t * sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd * delta_t;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a*delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd*delta_t;

		//write predicted sigma point into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}

	//print result
	//std::cout << "Xsig_pred = " << std::endl << Xsig_pred_ << std::endl;

	/*******************************************************************************
	* Predict mean and covariance
	******************************************************************************/
	//Lesson 7, section 24

	// set weights
	double weight_0 = lambda / (lambda + n_aug_);
	weights_(0) = weight_0;
	for (int i = 1; i<2 * n_aug_ + 1; i++) {  //2n+1 weights
		double weight = 0.5 / (n_aug_ + lambda);
		weights_(i) = weight;
	}

	//predicted state mean
	x_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}

	//std::cout << "predicted state out by bakr: " << std::endl;

	//predicted state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
											   // state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		//cout << "inside loop belal1\n";
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		//cout << "inside loop belal2\n";

		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;
		//cout << "inside loop belal3\n";

		P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
	}

	//print result
	/*std::cout << "Predicted state" << std::endl;
	std::cout << x_ << std::endl;
	std::cout << "Predicted covariance matrix" << std::endl;
	std::cout << P_ << std::endl;*/
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  //measurement matrix
	MatrixXd H_laser_ = MatrixXd(2, 5);
	H_laser_ << 1, 0, 0, 0, 0,
				0, 1, 0, 0, 0;

	//measurement covariance matrix - laser
	//This matrix must be given by the manufacturer
	MatrixXd R_laser_ = MatrixXd(2, 2);
	R_laser_ << std_laspx_, 0,
				0, std_laspx_;

	VectorXd z = meas_package.raw_measurements_;
	VectorXd z_pred = H_laser_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_laser_.transpose();
	MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	VectorXd new_x = x_ + (K * y);
	if (isReasonableNewX(new_x)) {
		x_ = new_x;
		long x_size = x_.size();
		MatrixXd I = MatrixXd::Identity(x_size, x_size);
		P_ = (I - K * H_laser_) * P_;
	}
}

bool UKF::isReasonableNewX(const VectorXd & newX)
{
	static int number_of_readings = 0;
	//ignore the first 4 calls to this function
	if (number_of_readings < 4)
	{
		number_of_readings++;
		return true;
	}
	const double new_x_val = newX(0);
	const double new_y_val = newX(1);
	const double new_xvel_val = newX(2);
	const double new_yvel_val = newX(3);

	const double curr_x_val = x_(0);
	const double curr_y_val = x_(1);
	const double curr_xvel_val = x_(2);
	const double curr_yvel_val = x_(3);

	if (abs(new_x_val - curr_x_val) > 0.3)return false;
	if (abs(new_y_val - curr_y_val) > 0.3)return false;
	/*if (abs(new_xvel_val - curr_xvel_val) > 0.6)return false;
	if (abs(new_yvel_val - curr_yvel_val) > 0.6)return false;*/

	return true;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  /*******************************************************************************
  * Predict Radar Sigma Points
  ******************************************************************************/
	//Lesson 7, section 27.

  //set measurement dimension, radar can measure r, phi, and r_dot
	int n_z = 3;

	//define spreading parameter
	double lambda = 3 - n_aug_;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

											   // extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                        //r
		Zsig(1, i) = atan2(p_y, p_x);                                 //phi
		Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//innovation covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
											   //residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_radr_*std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;
	S = S + R;

	//print result
	//std::cout << "z_pred: " << std::endl << z_pred << std::endl;
	//std::cout << "S: " << std::endl << S << std::endl;


	/*******************************************************************************
	* Update Radar
	******************************************************************************/
	//Lesson 7, section 30.

	VectorXd z = meas_package.raw_measurements_;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	/*******************************************************************************
	* Student part begin
	******************************************************************************/

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

											   //residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}

	//Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//residual
	VectorXd z_diff = z - z_pred;

	//angle normalization
	while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

	//update state mean and covariance matrix
	VectorXd newX = x_ + K * z_diff;
	if (isReasonableNewX(newX)) {
		x_ = newX;
	    P_ = P_ - K*S*K.transpose();
	}
	/*******************************************************************************
	* Student part end
	******************************************************************************/

	//print result
	/*std::cout << "Updated state x: " << std::endl << x_ << std::endl;
	std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;*/
}