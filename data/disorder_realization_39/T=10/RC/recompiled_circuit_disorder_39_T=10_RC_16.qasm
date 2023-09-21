OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0948148) q[0];
sx q[0];
rz(4.2098213) q[0];
sx q[0];
rz(9.8888483) q[0];
rz(-1.1820473) q[1];
sx q[1];
rz(-3.074475) q[1];
sx q[1];
rz(1.8571412) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0604295) q[0];
sx q[0];
rz(-1.2788749) q[0];
sx q[0];
rz(1.97776) q[0];
rz(-3.0797144) q[2];
sx q[2];
rz(-0.49052325) q[2];
sx q[2];
rz(-0.86531901) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7317036) q[1];
sx q[1];
rz(-1.1176425) q[1];
sx q[1];
rz(-1.1660006) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.554261) q[3];
sx q[3];
rz(-1.7214805) q[3];
sx q[3];
rz(2.2076803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7444732) q[2];
sx q[2];
rz(-2.2019272) q[2];
sx q[2];
rz(-2.822067) q[2];
rz(-0.5784936) q[3];
sx q[3];
rz(-2.6632023) q[3];
sx q[3];
rz(-0.67392504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4085061) q[0];
sx q[0];
rz(-1.5717614) q[0];
sx q[0];
rz(-0.52655667) q[0];
rz(0.56354848) q[1];
sx q[1];
rz(-1.6105517) q[1];
sx q[1];
rz(-0.79663509) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5986886) q[0];
sx q[0];
rz(-1.2447378) q[0];
sx q[0];
rz(-0.90755264) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57467069) q[2];
sx q[2];
rz(-1.4029014) q[2];
sx q[2];
rz(1.6811973) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.14256515) q[1];
sx q[1];
rz(-1.5579281) q[1];
sx q[1];
rz(-0.60728118) q[1];
rz(-pi) q[2];
rz(0.91005743) q[3];
sx q[3];
rz(-0.72477341) q[3];
sx q[3];
rz(-1.7006601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9600296) q[2];
sx q[2];
rz(-1.7627565) q[2];
sx q[2];
rz(-2.3550418) q[2];
rz(-0.49318796) q[3];
sx q[3];
rz(-1.1816053) q[3];
sx q[3];
rz(0.44737679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2717473) q[0];
sx q[0];
rz(-0.87690502) q[0];
sx q[0];
rz(1.440381) q[0];
rz(-2.4213743) q[1];
sx q[1];
rz(-0.59599382) q[1];
sx q[1];
rz(0.70297855) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7754606) q[0];
sx q[0];
rz(-1.4538987) q[0];
sx q[0];
rz(-1.2890105) q[0];
rz(1.9269283) q[2];
sx q[2];
rz(-2.0424358) q[2];
sx q[2];
rz(0.65442649) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30333334) q[1];
sx q[1];
rz(-1.6175744) q[1];
sx q[1];
rz(0.10951885) q[1];
rz(-pi) q[2];
rz(-0.17351563) q[3];
sx q[3];
rz(-2.1870038) q[3];
sx q[3];
rz(-0.71615744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8748223) q[2];
sx q[2];
rz(-1.4845994) q[2];
sx q[2];
rz(2.2300569) q[2];
rz(2.2276145) q[3];
sx q[3];
rz(-1.711288) q[3];
sx q[3];
rz(0.82733697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5660969) q[0];
sx q[0];
rz(-0.63593447) q[0];
sx q[0];
rz(0.77600586) q[0];
rz(1.2674747) q[1];
sx q[1];
rz(-2.0327366) q[1];
sx q[1];
rz(-2.5783096) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.027772) q[0];
sx q[0];
rz(-1.2978683) q[0];
sx q[0];
rz(0.91822894) q[0];
x q[1];
rz(-1.9070542) q[2];
sx q[2];
rz(-0.90869892) q[2];
sx q[2];
rz(1.3627571) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1174417) q[1];
sx q[1];
rz(-1.5772595) q[1];
sx q[1];
rz(-2.8698189) q[1];
rz(-pi) q[2];
rz(-1.8654278) q[3];
sx q[3];
rz(-2.2664824) q[3];
sx q[3];
rz(2.5254315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.48878601) q[2];
sx q[2];
rz(-2.2061009) q[2];
sx q[2];
rz(-2.9768067) q[2];
rz(0.22848836) q[3];
sx q[3];
rz(-0.27992862) q[3];
sx q[3];
rz(-0.94648186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8673458) q[0];
sx q[0];
rz(-2.2773401) q[0];
sx q[0];
rz(2.2891323) q[0];
rz(-2.7903941) q[1];
sx q[1];
rz(-0.57792592) q[1];
sx q[1];
rz(-1.9794827) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6974555) q[0];
sx q[0];
rz(-1.9636969) q[0];
sx q[0];
rz(-0.59209728) q[0];
rz(-2.2010872) q[2];
sx q[2];
rz(-0.61912196) q[2];
sx q[2];
rz(-1.3319912) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.64486849) q[1];
sx q[1];
rz(-2.2779896) q[1];
sx q[1];
rz(-2.1281388) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1245072) q[3];
sx q[3];
rz(-2.4199329) q[3];
sx q[3];
rz(0.34860308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6953485) q[2];
sx q[2];
rz(-2.7480795) q[2];
sx q[2];
rz(-1.8943141) q[2];
rz(-0.013109664) q[3];
sx q[3];
rz(-1.7316827) q[3];
sx q[3];
rz(-2.9124027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2728249) q[0];
sx q[0];
rz(-1.8936963) q[0];
sx q[0];
rz(-0.75063467) q[0];
rz(1.1095095) q[1];
sx q[1];
rz(-0.37056071) q[1];
sx q[1];
rz(-1.8355339) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36201492) q[0];
sx q[0];
rz(-2.0834196) q[0];
sx q[0];
rz(-2.6984452) q[0];
rz(2.7517031) q[2];
sx q[2];
rz(-0.49553686) q[2];
sx q[2];
rz(-1.3930266) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8667824) q[1];
sx q[1];
rz(-0.85046235) q[1];
sx q[1];
rz(-1.5809098) q[1];
rz(0.2209729) q[3];
sx q[3];
rz(-0.60022012) q[3];
sx q[3];
rz(-2.5252987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3874454) q[2];
sx q[2];
rz(-1.0070609) q[2];
sx q[2];
rz(0.60069096) q[2];
rz(-1.0026275) q[3];
sx q[3];
rz(-2.6032344) q[3];
sx q[3];
rz(2.7887662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7464741) q[0];
sx q[0];
rz(-1.6858608) q[0];
sx q[0];
rz(-1.0466928) q[0];
rz(-1.612161) q[1];
sx q[1];
rz(-1.7144831) q[1];
sx q[1];
rz(-2.7244862) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1313275) q[0];
sx q[0];
rz(-1.7908887) q[0];
sx q[0];
rz(3.1165645) q[0];
rz(-pi) q[1];
rz(-0.30714005) q[2];
sx q[2];
rz(-1.2899613) q[2];
sx q[2];
rz(1.839523) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.68361359) q[1];
sx q[1];
rz(-2.2375467) q[1];
sx q[1];
rz(-0.97138202) q[1];
x q[2];
rz(-2.4884175) q[3];
sx q[3];
rz(-1.8110868) q[3];
sx q[3];
rz(-2.0865292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.0059011857) q[2];
sx q[2];
rz(-0.4726755) q[2];
sx q[2];
rz(-1.3674412) q[2];
rz(-0.55316365) q[3];
sx q[3];
rz(-1.4033214) q[3];
sx q[3];
rz(0.52545351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.816514) q[0];
sx q[0];
rz(-1.8186318) q[0];
sx q[0];
rz(-1.9518071) q[0];
rz(-1.4272383) q[1];
sx q[1];
rz(-0.27806565) q[1];
sx q[1];
rz(0.11238012) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6067137) q[0];
sx q[0];
rz(-2.0438072) q[0];
sx q[0];
rz(-2.9472449) q[0];
x q[1];
rz(-0.24765315) q[2];
sx q[2];
rz(-2.7981749) q[2];
sx q[2];
rz(-2.0695956) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1689414) q[1];
sx q[1];
rz(-2.2691138) q[1];
sx q[1];
rz(-2.6886743) q[1];
rz(-pi) q[2];
rz(1.3120679) q[3];
sx q[3];
rz(-1.6885763) q[3];
sx q[3];
rz(-0.044101276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0671976) q[2];
sx q[2];
rz(-1.9125568) q[2];
sx q[2];
rz(-1.3486264) q[2];
rz(-1.9366692) q[3];
sx q[3];
rz(-1.4068853) q[3];
sx q[3];
rz(2.9437734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39712054) q[0];
sx q[0];
rz(-1.67698) q[0];
sx q[0];
rz(-0.034974139) q[0];
rz(2.2947625) q[1];
sx q[1];
rz(-1.7085107) q[1];
sx q[1];
rz(2.2299178) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23110403) q[0];
sx q[0];
rz(-1.3470874) q[0];
sx q[0];
rz(1.847812) q[0];
x q[1];
rz(-0.54347221) q[2];
sx q[2];
rz(-1.4620355) q[2];
sx q[2];
rz(2.6639338) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1940143) q[1];
sx q[1];
rz(-1.7742426) q[1];
sx q[1];
rz(-0.46781637) q[1];
rz(-pi) q[2];
rz(-1.0802286) q[3];
sx q[3];
rz(-2.053223) q[3];
sx q[3];
rz(-0.8182984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4336865) q[2];
sx q[2];
rz(-0.5138548) q[2];
sx q[2];
rz(-0.24924499) q[2];
rz(0.76672673) q[3];
sx q[3];
rz(-1.8079575) q[3];
sx q[3];
rz(-2.7419817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7837759) q[0];
sx q[0];
rz(-2.0932842) q[0];
sx q[0];
rz(-2.7695079) q[0];
rz(2.5601939) q[1];
sx q[1];
rz(-2.364368) q[1];
sx q[1];
rz(1.4153597) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6045195) q[0];
sx q[0];
rz(-1.8329289) q[0];
sx q[0];
rz(1.7781236) q[0];
rz(2.084311) q[2];
sx q[2];
rz(-1.0143177) q[2];
sx q[2];
rz(-1.8884115) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.806554) q[1];
sx q[1];
rz(-1.6400669) q[1];
sx q[1];
rz(-3.0663504) q[1];
rz(-0.98202242) q[3];
sx q[3];
rz(-2.5128799) q[3];
sx q[3];
rz(-0.94129896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.32594484) q[2];
sx q[2];
rz(-0.14020136) q[2];
sx q[2];
rz(-2.9369205) q[2];
rz(-1.7278016) q[3];
sx q[3];
rz(-1.9982343) q[3];
sx q[3];
rz(2.0457101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6407912) q[0];
sx q[0];
rz(-1.6727819) q[0];
sx q[0];
rz(0.93760136) q[0];
rz(-1.5851371) q[1];
sx q[1];
rz(-1.382348) q[1];
sx q[1];
rz(1.4437645) q[1];
rz(1.1164222) q[2];
sx q[2];
rz(-1.6797671) q[2];
sx q[2];
rz(-1.4327008) q[2];
rz(0.99863573) q[3];
sx q[3];
rz(-1.640366) q[3];
sx q[3];
rz(0.58983005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];