OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62701464) q[0];
sx q[0];
rz(3.7778683) q[0];
sx q[0];
rz(10.502622) q[0];
rz(1.0765422) q[1];
sx q[1];
rz(-0.62907469) q[1];
sx q[1];
rz(-1.0318626) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54993486) q[0];
sx q[0];
rz(-3.1374133) q[0];
sx q[0];
rz(-3.0463534) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3276153) q[2];
sx q[2];
rz(-1.2300101) q[2];
sx q[2];
rz(1.5208517) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.061432211) q[1];
sx q[1];
rz(-2.5974063) q[1];
sx q[1];
rz(-2.4341492) q[1];
x q[2];
rz(0.21764619) q[3];
sx q[3];
rz(-0.34322327) q[3];
sx q[3];
rz(1.1962593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91118139) q[2];
sx q[2];
rz(-1.958467) q[2];
sx q[2];
rz(-2.6739056) q[2];
rz(2.6905401) q[3];
sx q[3];
rz(-0.39877287) q[3];
sx q[3];
rz(0.27802813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4462047) q[0];
sx q[0];
rz(-2.9187293) q[0];
sx q[0];
rz(-0.15431246) q[0];
rz(2.800324) q[1];
sx q[1];
rz(-2.7879265) q[1];
sx q[1];
rz(2.7214859) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4898191) q[0];
sx q[0];
rz(-0.66733995) q[0];
sx q[0];
rz(0.31489189) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2809199) q[2];
sx q[2];
rz(-1.6063074) q[2];
sx q[2];
rz(1.4461609) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5771835) q[1];
sx q[1];
rz(-2.0415876) q[1];
sx q[1];
rz(-1.5723349) q[1];
rz(2.9734334) q[3];
sx q[3];
rz(-0.75440591) q[3];
sx q[3];
rz(-1.8969632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5306065) q[2];
sx q[2];
rz(-2.2425118) q[2];
sx q[2];
rz(2.3685624) q[2];
rz(-2.2981339) q[3];
sx q[3];
rz(-1.5777028) q[3];
sx q[3];
rz(2.5696866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.968349) q[0];
sx q[0];
rz(-1.0272212) q[0];
sx q[0];
rz(-2.9395043) q[0];
rz(-0.48187065) q[1];
sx q[1];
rz(-0.20129573) q[1];
sx q[1];
rz(-2.2352107) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2807686) q[0];
sx q[0];
rz(-1.1204136) q[0];
sx q[0];
rz(0.76669873) q[0];
rz(-pi) q[1];
rz(2.2041476) q[2];
sx q[2];
rz(-0.54493947) q[2];
sx q[2];
rz(2.3505806) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6828595) q[1];
sx q[1];
rz(-0.56536973) q[1];
sx q[1];
rz(-0.302178) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0624978) q[3];
sx q[3];
rz(-2.8064499) q[3];
sx q[3];
rz(2.3329874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.071094461) q[2];
sx q[2];
rz(-1.1134032) q[2];
sx q[2];
rz(0.58007288) q[2];
rz(-1.5602559) q[3];
sx q[3];
rz(-1.7507078) q[3];
sx q[3];
rz(1.7177379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4105014) q[0];
sx q[0];
rz(-3.0803362) q[0];
sx q[0];
rz(0.047792338) q[0];
rz(1.8675249) q[1];
sx q[1];
rz(-1.85227) q[1];
sx q[1];
rz(3.0963669) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84789373) q[0];
sx q[0];
rz(-1.4220823) q[0];
sx q[0];
rz(2.2151164) q[0];
rz(-pi) q[1];
rz(0.93542807) q[2];
sx q[2];
rz(-1.3029939) q[2];
sx q[2];
rz(-1.2768465) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9103376) q[1];
sx q[1];
rz(-1.567776) q[1];
sx q[1];
rz(-1.5652204) q[1];
rz(-pi) q[2];
rz(0.91462709) q[3];
sx q[3];
rz(-2.4883929) q[3];
sx q[3];
rz(0.40756215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4101326) q[2];
sx q[2];
rz(-1.2612017) q[2];
sx q[2];
rz(0.88978466) q[2];
rz(2.5951923) q[3];
sx q[3];
rz(-1.0011287) q[3];
sx q[3];
rz(-2.8304097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6379717) q[0];
sx q[0];
rz(-1.4542955) q[0];
sx q[0];
rz(1.7244435) q[0];
rz(-1.1196989) q[1];
sx q[1];
rz(-1.3816625) q[1];
sx q[1];
rz(-1.959257) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6851824) q[0];
sx q[0];
rz(-1.6401924) q[0];
sx q[0];
rz(0.38763898) q[0];
rz(-pi) q[1];
rz(-2.3619217) q[2];
sx q[2];
rz(-0.60616772) q[2];
sx q[2];
rz(0.066628284) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.87403832) q[1];
sx q[1];
rz(-2.2910765) q[1];
sx q[1];
rz(-1.6362067) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0145217) q[3];
sx q[3];
rz(-2.2699353) q[3];
sx q[3];
rz(-2.3026031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.22415367) q[2];
sx q[2];
rz(-2.3771693) q[2];
sx q[2];
rz(2.7177641) q[2];
rz(2.0297) q[3];
sx q[3];
rz(-0.49911505) q[3];
sx q[3];
rz(-0.65139884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79913419) q[0];
sx q[0];
rz(-0.0722216) q[0];
sx q[0];
rz(-2.1891731) q[0];
rz(2.3656288) q[1];
sx q[1];
rz(-0.9351848) q[1];
sx q[1];
rz(-1.5663358) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.562378) q[0];
sx q[0];
rz(-1.5327306) q[0];
sx q[0];
rz(-1.5910351) q[0];
rz(1.6265154) q[2];
sx q[2];
rz(-2.5341659) q[2];
sx q[2];
rz(-0.78616963) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.414622) q[1];
sx q[1];
rz(-1.7900253) q[1];
sx q[1];
rz(0.61466932) q[1];
rz(2.0616343) q[3];
sx q[3];
rz(-1.4075402) q[3];
sx q[3];
rz(-1.8051683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39885193) q[2];
sx q[2];
rz(-0.80790085) q[2];
sx q[2];
rz(-2.3952132) q[2];
rz(2.0505203) q[3];
sx q[3];
rz(-1.4336136) q[3];
sx q[3];
rz(2.5892042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2090476) q[0];
sx q[0];
rz(-0.046367558) q[0];
sx q[0];
rz(2.5982502) q[0];
rz(0.53687334) q[1];
sx q[1];
rz(-0.32506341) q[1];
sx q[1];
rz(0.12319014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9730102) q[0];
sx q[0];
rz(-1.5401398) q[0];
sx q[0];
rz(-0.93498303) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4692612) q[2];
sx q[2];
rz(-0.62158442) q[2];
sx q[2];
rz(0.42447916) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0052332) q[1];
sx q[1];
rz(-0.87456223) q[1];
sx q[1];
rz(-0.74058786) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76456548) q[3];
sx q[3];
rz(-0.95029059) q[3];
sx q[3];
rz(-2.4303049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8366375) q[2];
sx q[2];
rz(-1.7540437) q[2];
sx q[2];
rz(-0.48509625) q[2];
rz(-1.9890316) q[3];
sx q[3];
rz(-1.3644812) q[3];
sx q[3];
rz(-0.047957234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7570067) q[0];
sx q[0];
rz(-2.204019) q[0];
sx q[0];
rz(-2.6458929) q[0];
rz(-0.063830201) q[1];
sx q[1];
rz(-0.7494691) q[1];
sx q[1];
rz(-3.0705423) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97489965) q[0];
sx q[0];
rz(-1.3314795) q[0];
sx q[0];
rz(-0.47310754) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60511968) q[2];
sx q[2];
rz(-1.7150884) q[2];
sx q[2];
rz(-0.52478204) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7361765) q[1];
sx q[1];
rz(-1.725493) q[1];
sx q[1];
rz(-1.8150034) q[1];
rz(2.9661353) q[3];
sx q[3];
rz(-2.6883295) q[3];
sx q[3];
rz(-2.3153265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0157328) q[2];
sx q[2];
rz(-1.6682397) q[2];
sx q[2];
rz(-2.2802172) q[2];
rz(1.8371948) q[3];
sx q[3];
rz(-2.5671037) q[3];
sx q[3];
rz(1.5484352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9651589) q[0];
sx q[0];
rz(-0.49089828) q[0];
sx q[0];
rz(1.4392256) q[0];
rz(0.08055117) q[1];
sx q[1];
rz(-1.4713902) q[1];
sx q[1];
rz(1.390994) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8497365) q[0];
sx q[0];
rz(-2.448659) q[0];
sx q[0];
rz(-0.87397184) q[0];
rz(-1.2660145) q[2];
sx q[2];
rz(-1.5902747) q[2];
sx q[2];
rz(1.2661042) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1999368) q[1];
sx q[1];
rz(-2.1547531) q[1];
sx q[1];
rz(0.55364174) q[1];
x q[2];
rz(-1.9506878) q[3];
sx q[3];
rz(-1.2244488) q[3];
sx q[3];
rz(-2.430918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78802687) q[2];
sx q[2];
rz(-1.8466419) q[2];
sx q[2];
rz(-3.0143152) q[2];
rz(0.92653972) q[3];
sx q[3];
rz(-0.24181952) q[3];
sx q[3];
rz(1.4827137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7552898) q[0];
sx q[0];
rz(-0.64978623) q[0];
sx q[0];
rz(1.5007098) q[0];
rz(0.81037784) q[1];
sx q[1];
rz(-2.2893298) q[1];
sx q[1];
rz(-2.3694029) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92896748) q[0];
sx q[0];
rz(-1.0662931) q[0];
sx q[0];
rz(2.3113768) q[0];
rz(-pi) q[1];
rz(0.61913808) q[2];
sx q[2];
rz(-0.35751128) q[2];
sx q[2];
rz(-0.50019568) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4160847) q[1];
sx q[1];
rz(-0.5495175) q[1];
sx q[1];
rz(-1.548442) q[1];
x q[2];
rz(1.6343035) q[3];
sx q[3];
rz(-2.4606703) q[3];
sx q[3];
rz(-1.6977967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3810252) q[2];
sx q[2];
rz(-2.3367391) q[2];
sx q[2];
rz(0.62925657) q[2];
rz(0.73090807) q[3];
sx q[3];
rz(-2.2980502) q[3];
sx q[3];
rz(-1.4430911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6998941) q[0];
sx q[0];
rz(-2.6581673) q[0];
sx q[0];
rz(2.7375258) q[0];
rz(-0.40456698) q[1];
sx q[1];
rz(-1.7351983) q[1];
sx q[1];
rz(2.4529967) q[1];
rz(2.9108638) q[2];
sx q[2];
rz(-0.28906549) q[2];
sx q[2];
rz(-1.1171823) q[2];
rz(2.6062905) q[3];
sx q[3];
rz(-1.3131159) q[3];
sx q[3];
rz(2.5406607) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
