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
rz(3.0313015) q[0];
sx q[0];
rz(-1.8561441) q[0];
sx q[0];
rz(0.31874803) q[0];
rz(4.3238001) q[1];
sx q[1];
rz(4.8053513) q[1];
sx q[1];
rz(4.3045192) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010745239) q[0];
sx q[0];
rz(-1.3358572) q[0];
sx q[0];
rz(1.9128591) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0672795) q[2];
sx q[2];
rz(-0.19467672) q[2];
sx q[2];
rz(-0.014460221) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0981465) q[1];
sx q[1];
rz(-0.43546989) q[1];
sx q[1];
rz(0.79186073) q[1];
rz(-0.57489328) q[3];
sx q[3];
rz(-1.5539317) q[3];
sx q[3];
rz(-0.3080388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17243871) q[2];
sx q[2];
rz(-2.6814851) q[2];
sx q[2];
rz(-3.0058506) q[2];
rz(0.33484778) q[3];
sx q[3];
rz(-1.2232774) q[3];
sx q[3];
rz(-0.39723435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44884509) q[0];
sx q[0];
rz(-1.0053585) q[0];
sx q[0];
rz(2.1951065) q[0];
rz(1.954156) q[1];
sx q[1];
rz(-0.58381909) q[1];
sx q[1];
rz(2.8576287) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07325726) q[0];
sx q[0];
rz(-0.26480246) q[0];
sx q[0];
rz(-2.3786484) q[0];
rz(-0.87209629) q[2];
sx q[2];
rz(-1.3882033) q[2];
sx q[2];
rz(-2.3645949) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3888106) q[1];
sx q[1];
rz(-1.7080664) q[1];
sx q[1];
rz(0.42690446) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3881298) q[3];
sx q[3];
rz(-1.2468657) q[3];
sx q[3];
rz(0.38979724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2076608) q[2];
sx q[2];
rz(-1.8457103) q[2];
sx q[2];
rz(2.3409823) q[2];
rz(-1.8396395) q[3];
sx q[3];
rz(-2.0079565) q[3];
sx q[3];
rz(0.058066644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0170853) q[0];
sx q[0];
rz(-0.43211102) q[0];
sx q[0];
rz(2.3486163) q[0];
rz(-0.68603459) q[1];
sx q[1];
rz(-1.7162836) q[1];
sx q[1];
rz(-2.3763903) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52872872) q[0];
sx q[0];
rz(-2.9224797) q[0];
sx q[0];
rz(-1.8414453) q[0];
rz(0.051570895) q[2];
sx q[2];
rz(-1.602939) q[2];
sx q[2];
rz(-0.78173897) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0211693) q[1];
sx q[1];
rz(-0.83008728) q[1];
sx q[1];
rz(2.7629025) q[1];
rz(2.5038165) q[3];
sx q[3];
rz(-2.2368663) q[3];
sx q[3];
rz(3.1036669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5746295) q[2];
sx q[2];
rz(-0.69090635) q[2];
sx q[2];
rz(-1.2464657) q[2];
rz(-1.9555107) q[3];
sx q[3];
rz(-1.7639152) q[3];
sx q[3];
rz(2.932909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3935811) q[0];
sx q[0];
rz(-1.8210541) q[0];
sx q[0];
rz(3.0082974) q[0];
rz(-2.8721299) q[1];
sx q[1];
rz(-0.91454426) q[1];
sx q[1];
rz(2.6752313) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8036706) q[0];
sx q[0];
rz(-0.9804157) q[0];
sx q[0];
rz(-1.7147816) q[0];
rz(-0.94302098) q[2];
sx q[2];
rz(-1.7258712) q[2];
sx q[2];
rz(2.6076406) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.129648) q[1];
sx q[1];
rz(-2.5135871) q[1];
sx q[1];
rz(0.25917128) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5047938) q[3];
sx q[3];
rz(-2.0836692) q[3];
sx q[3];
rz(-2.6309225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2103297) q[2];
sx q[2];
rz(-0.42079058) q[2];
sx q[2];
rz(-2.8974864) q[2];
rz(-1.7804451) q[3];
sx q[3];
rz(-2.1972392) q[3];
sx q[3];
rz(-1.2066427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72531438) q[0];
sx q[0];
rz(-0.84753528) q[0];
sx q[0];
rz(3.0552979) q[0];
rz(-1.7591954) q[1];
sx q[1];
rz(-2.081213) q[1];
sx q[1];
rz(1.28654) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79199687) q[0];
sx q[0];
rz(-1.9583869) q[0];
sx q[0];
rz(0.18328615) q[0];
rz(-pi) q[1];
rz(2.9164739) q[2];
sx q[2];
rz(-0.79155603) q[2];
sx q[2];
rz(-2.5560372) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4205192) q[1];
sx q[1];
rz(-0.87158485) q[1];
sx q[1];
rz(-3.0341604) q[1];
rz(-pi) q[2];
rz(2.6161353) q[3];
sx q[3];
rz(-0.94949978) q[3];
sx q[3];
rz(2.3182825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4061819) q[2];
sx q[2];
rz(-1.7291131) q[2];
sx q[2];
rz(2.0162876) q[2];
rz(2.3451037) q[3];
sx q[3];
rz(-0.73553604) q[3];
sx q[3];
rz(-0.19851941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7859802) q[0];
sx q[0];
rz(-2.8481843) q[0];
sx q[0];
rz(-2.8134213) q[0];
rz(-2.7525821) q[1];
sx q[1];
rz(-1.5711454) q[1];
sx q[1];
rz(-1.4555812) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4220912) q[0];
sx q[0];
rz(-1.6515284) q[0];
sx q[0];
rz(-0.21105612) q[0];
rz(-0.12673086) q[2];
sx q[2];
rz(-1.4369082) q[2];
sx q[2];
rz(2.9478879) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.38992369) q[1];
sx q[1];
rz(-1.3397386) q[1];
sx q[1];
rz(2.2242145) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1001411) q[3];
sx q[3];
rz(-1.2123479) q[3];
sx q[3];
rz(0.35669092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1652611) q[2];
sx q[2];
rz(-0.92504048) q[2];
sx q[2];
rz(-2.6791005) q[2];
rz(-1.0235323) q[3];
sx q[3];
rz(-2.0868802) q[3];
sx q[3];
rz(-0.83824497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3384712) q[0];
sx q[0];
rz(-1.4305038) q[0];
sx q[0];
rz(-2.9679003) q[0];
rz(0.22131418) q[1];
sx q[1];
rz(-0.68562713) q[1];
sx q[1];
rz(-1.7783222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1376138) q[0];
sx q[0];
rz(-2.6948746) q[0];
sx q[0];
rz(2.0017308) q[0];
rz(-pi) q[1];
rz(-1.4837711) q[2];
sx q[2];
rz(-1.489893) q[2];
sx q[2];
rz(1.9593585) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.89243556) q[1];
sx q[1];
rz(-1.1022864) q[1];
sx q[1];
rz(-2.4484642) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24226455) q[3];
sx q[3];
rz(-1.8349832) q[3];
sx q[3];
rz(1.3380877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26019105) q[2];
sx q[2];
rz(-1.1611725) q[2];
sx q[2];
rz(1.6960404) q[2];
rz(2.3675303) q[3];
sx q[3];
rz(-2.4249707) q[3];
sx q[3];
rz(-1.6707481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.562302) q[0];
sx q[0];
rz(-0.93920541) q[0];
sx q[0];
rz(-0.24765177) q[0];
rz(0.66894764) q[1];
sx q[1];
rz(-1.9525783) q[1];
sx q[1];
rz(0.98562366) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0730608) q[0];
sx q[0];
rz(-1.685623) q[0];
sx q[0];
rz(-0.5590183) q[0];
rz(-pi) q[1];
rz(2.5778887) q[2];
sx q[2];
rz(-1.3022175) q[2];
sx q[2];
rz(-0.33252663) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0874153) q[1];
sx q[1];
rz(-1.8761411) q[1];
sx q[1];
rz(1.9806238) q[1];
rz(-0.0019916742) q[3];
sx q[3];
rz(-2.5466306) q[3];
sx q[3];
rz(1.4548986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.140427) q[2];
sx q[2];
rz(-2.3900034) q[2];
sx q[2];
rz(1.8482194) q[2];
rz(1.2398531) q[3];
sx q[3];
rz(-0.69609061) q[3];
sx q[3];
rz(0.70824879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2769315) q[0];
sx q[0];
rz(-1.7129352) q[0];
sx q[0];
rz(-0.96555936) q[0];
rz(-2.7652265) q[1];
sx q[1];
rz(-1.3220359) q[1];
sx q[1];
rz(-1.2300864) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.108889) q[0];
sx q[0];
rz(-2.8710693) q[0];
sx q[0];
rz(-2.9530086) q[0];
x q[1];
rz(-1.240991) q[2];
sx q[2];
rz(-1.9422741) q[2];
sx q[2];
rz(-0.35201752) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2918307) q[1];
sx q[1];
rz(-0.9660143) q[1];
sx q[1];
rz(-1.6563708) q[1];
x q[2];
rz(2.3283703) q[3];
sx q[3];
rz(-1.5322161) q[3];
sx q[3];
rz(0.67594066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2824715) q[2];
sx q[2];
rz(-2.4716447) q[2];
sx q[2];
rz(0.13151375) q[2];
rz(-0.48111835) q[3];
sx q[3];
rz(-0.93004623) q[3];
sx q[3];
rz(0.49829811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1199101) q[0];
sx q[0];
rz(-2.5639738) q[0];
sx q[0];
rz(-2.6525894) q[0];
rz(-1.3267714) q[1];
sx q[1];
rz(-1.2811456) q[1];
sx q[1];
rz(-0.16924032) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3796462) q[0];
sx q[0];
rz(-0.45416203) q[0];
sx q[0];
rz(1.2962925) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6851114) q[2];
sx q[2];
rz(-0.20083961) q[2];
sx q[2];
rz(2.0690837) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9368912) q[1];
sx q[1];
rz(-1.1549779) q[1];
sx q[1];
rz(-2.4067626) q[1];
x q[2];
rz(1.2792688) q[3];
sx q[3];
rz(-1.3025373) q[3];
sx q[3];
rz(-0.30941468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5883611) q[2];
sx q[2];
rz(-2.0560122) q[2];
sx q[2];
rz(0.21151839) q[2];
rz(-1.616098) q[3];
sx q[3];
rz(-1.0001405) q[3];
sx q[3];
rz(0.26743993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78032988) q[0];
sx q[0];
rz(-1.4365256) q[0];
sx q[0];
rz(-2.8509675) q[0];
rz(0.81193874) q[1];
sx q[1];
rz(-0.62807905) q[1];
sx q[1];
rz(-1.5608578) q[1];
rz(-0.43368922) q[2];
sx q[2];
rz(-1.5530752) q[2];
sx q[2];
rz(-1.8032522) q[2];
rz(-2.5200882) q[3];
sx q[3];
rz(-2.2198698) q[3];
sx q[3];
rz(2.6563016) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
