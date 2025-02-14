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
rz(-2.8228446) q[0];
rz(1.1822074) q[1];
sx q[1];
rz(-1.6637586) q[1];
sx q[1];
rz(-1.1629265) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1308474) q[0];
sx q[0];
rz(-1.8057355) q[0];
sx q[0];
rz(1.9128591) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0672795) q[2];
sx q[2];
rz(-2.9469159) q[2];
sx q[2];
rz(-0.014460221) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.27027495) q[1];
sx q[1];
rz(-1.2658889) q[1];
sx q[1];
rz(-0.31590806) q[1];
x q[2];
rz(-1.5908904) q[3];
sx q[3];
rz(-2.1455975) q[3];
sx q[3];
rz(-1.8679096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9691539) q[2];
sx q[2];
rz(-0.46010751) q[2];
sx q[2];
rz(0.13574204) q[2];
rz(0.33484778) q[3];
sx q[3];
rz(-1.9183153) q[3];
sx q[3];
rz(-2.7443583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.44884509) q[0];
sx q[0];
rz(-2.1362342) q[0];
sx q[0];
rz(-0.94648615) q[0];
rz(-1.954156) q[1];
sx q[1];
rz(-2.5577736) q[1];
sx q[1];
rz(-0.28396398) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2427802) q[0];
sx q[0];
rz(-1.3889342) q[0];
sx q[0];
rz(-2.9480431) q[0];
rz(-pi) q[1];
rz(-0.23663972) q[2];
sx q[2];
rz(-2.2556117) q[2];
sx q[2];
rz(0.94517148) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.52580035) q[1];
sx q[1];
rz(-2.6944571) q[1];
sx q[1];
rz(2.8195803) q[1];
x q[2];
rz(1.2228187) q[3];
sx q[3];
rz(-1.2038411) q[3];
sx q[3];
rz(1.8311794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9339319) q[2];
sx q[2];
rz(-1.2958823) q[2];
sx q[2];
rz(-0.80061039) q[2];
rz(-1.8396395) q[3];
sx q[3];
rz(-2.0079565) q[3];
sx q[3];
rz(-3.083526) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0170853) q[0];
sx q[0];
rz(-2.7094816) q[0];
sx q[0];
rz(-2.3486163) q[0];
rz(2.4555581) q[1];
sx q[1];
rz(-1.7162836) q[1];
sx q[1];
rz(0.7652024) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8350462) q[0];
sx q[0];
rz(-1.6289428) q[0];
sx q[0];
rz(1.7821728) q[0];
x q[1];
rz(1.6029818) q[2];
sx q[2];
rz(-1.6223406) q[2];
sx q[2];
rz(-0.79071617) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0211693) q[1];
sx q[1];
rz(-2.3115054) q[1];
sx q[1];
rz(-2.7629025) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3451557) q[3];
sx q[3];
rz(-2.0579866) q[3];
sx q[3];
rz(2.0381442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.56696314) q[2];
sx q[2];
rz(-2.4506863) q[2];
sx q[2];
rz(-1.2464657) q[2];
rz(1.9555107) q[3];
sx q[3];
rz(-1.7639152) q[3];
sx q[3];
rz(0.20868364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3935811) q[0];
sx q[0];
rz(-1.3205386) q[0];
sx q[0];
rz(0.13329521) q[0];
rz(-2.8721299) q[1];
sx q[1];
rz(-2.2270484) q[1];
sx q[1];
rz(-2.6752313) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5488777) q[0];
sx q[0];
rz(-0.60564954) q[0];
sx q[0];
rz(0.2109419) q[0];
x q[1];
rz(1.830929) q[2];
sx q[2];
rz(-2.4974653) q[2];
sx q[2];
rz(0.82714236) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0119446) q[1];
sx q[1];
rz(-2.5135871) q[1];
sx q[1];
rz(-2.8824214) q[1];
rz(2.6277873) q[3];
sx q[3];
rz(-1.5132959) q[3];
sx q[3];
rz(-2.1138885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93126297) q[2];
sx q[2];
rz(-2.7208021) q[2];
sx q[2];
rz(2.8974864) q[2];
rz(-1.3611475) q[3];
sx q[3];
rz(-0.94435349) q[3];
sx q[3];
rz(1.9349499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4162783) q[0];
sx q[0];
rz(-2.2940574) q[0];
sx q[0];
rz(-0.086294802) q[0];
rz(1.3823973) q[1];
sx q[1];
rz(-1.0603797) q[1];
sx q[1];
rz(-1.28654) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3495958) q[0];
sx q[0];
rz(-1.9583869) q[0];
sx q[0];
rz(2.9583065) q[0];
rz(-pi) q[1];
rz(-2.9164739) q[2];
sx q[2];
rz(-0.79155603) q[2];
sx q[2];
rz(2.5560372) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5865422) q[1];
sx q[1];
rz(-2.4355531) q[1];
sx q[1];
rz(1.6976188) q[1];
rz(-2.1819918) q[3];
sx q[3];
rz(-2.3510076) q[3];
sx q[3];
rz(-0.039856002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7354108) q[2];
sx q[2];
rz(-1.7291131) q[2];
sx q[2];
rz(-2.0162876) q[2];
rz(2.3451037) q[3];
sx q[3];
rz(-2.4060566) q[3];
sx q[3];
rz(0.19851941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7859802) q[0];
sx q[0];
rz(-0.2934083) q[0];
sx q[0];
rz(0.32817131) q[0];
rz(-0.38901058) q[1];
sx q[1];
rz(-1.5704472) q[1];
sx q[1];
rz(1.6860115) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7195014) q[0];
sx q[0];
rz(-1.4900643) q[0];
sx q[0];
rz(-2.9305365) q[0];
rz(-pi) q[1];
rz(-1.7057538) q[2];
sx q[2];
rz(-1.696387) q[2];
sx q[2];
rz(1.3940982) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1342866) q[1];
sx q[1];
rz(-2.2040229) q[1];
sx q[1];
rz(0.28805201) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4606299) q[3];
sx q[3];
rz(-0.36073437) q[3];
sx q[3];
rz(-2.9025789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1652611) q[2];
sx q[2];
rz(-2.2165522) q[2];
sx q[2];
rz(2.6791005) q[2];
rz(2.1180604) q[3];
sx q[3];
rz(-2.0868802) q[3];
sx q[3];
rz(-0.83824497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3384712) q[0];
sx q[0];
rz(-1.7110889) q[0];
sx q[0];
rz(2.9679003) q[0];
rz(-0.22131418) q[1];
sx q[1];
rz(-2.4559655) q[1];
sx q[1];
rz(-1.7783222) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46745978) q[0];
sx q[0];
rz(-1.1674351) q[0];
sx q[0];
rz(2.9441071) q[0];
rz(-2.3214746) q[2];
sx q[2];
rz(-3.0228399) q[2];
sx q[2];
rz(1.135716) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9611836) q[1];
sx q[1];
rz(-0.81432691) q[1];
sx q[1];
rz(0.66988887) q[1];
x q[2];
rz(-2.8993281) q[3];
sx q[3];
rz(-1.8349832) q[3];
sx q[3];
rz(-1.3380877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8814016) q[2];
sx q[2];
rz(-1.9804201) q[2];
sx q[2];
rz(1.6960404) q[2];
rz(-2.3675303) q[3];
sx q[3];
rz(-0.71662199) q[3];
sx q[3];
rz(1.4708446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5792907) q[0];
sx q[0];
rz(-2.2023872) q[0];
sx q[0];
rz(2.8939409) q[0];
rz(-2.472645) q[1];
sx q[1];
rz(-1.9525783) q[1];
sx q[1];
rz(-2.155969) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0730608) q[0];
sx q[0];
rz(-1.685623) q[0];
sx q[0];
rz(-0.5590183) q[0];
rz(2.5778887) q[2];
sx q[2];
rz(-1.3022175) q[2];
sx q[2];
rz(2.809066) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0528763) q[1];
sx q[1];
rz(-0.50584882) q[1];
sx q[1];
rz(-2.2400676) q[1];
rz(-pi) q[2];
rz(1.5721442) q[3];
sx q[3];
rz(-0.97583563) q[3];
sx q[3];
rz(1.6890989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0011657) q[2];
sx q[2];
rz(-2.3900034) q[2];
sx q[2];
rz(-1.8482194) q[2];
rz(-1.2398531) q[3];
sx q[3];
rz(-2.445502) q[3];
sx q[3];
rz(-2.4333439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86466113) q[0];
sx q[0];
rz(-1.4286574) q[0];
sx q[0];
rz(0.96555936) q[0];
rz(0.37636617) q[1];
sx q[1];
rz(-1.8195567) q[1];
sx q[1];
rz(-1.9115062) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2282283) q[0];
sx q[0];
rz(-1.8364062) q[0];
sx q[0];
rz(1.5188541) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69355857) q[2];
sx q[2];
rz(-0.4916113) q[2];
sx q[2];
rz(-0.40406057) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.849762) q[1];
sx q[1];
rz(-2.1755784) q[1];
sx q[1];
rz(-1.4852219) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8132223) q[3];
sx q[3];
rz(-1.6093765) q[3];
sx q[3];
rz(-2.465652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2824715) q[2];
sx q[2];
rz(-0.66994795) q[2];
sx q[2];
rz(-3.0100789) q[2];
rz(2.6604743) q[3];
sx q[3];
rz(-2.2115464) q[3];
sx q[3];
rz(-0.49829811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1199101) q[0];
sx q[0];
rz(-0.57761884) q[0];
sx q[0];
rz(0.48900327) q[0];
rz(1.8148212) q[1];
sx q[1];
rz(-1.860447) q[1];
sx q[1];
rz(-2.9723523) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76194644) q[0];
sx q[0];
rz(-0.45416203) q[0];
sx q[0];
rz(1.8453002) q[0];
x q[1];
rz(1.4812977) q[2];
sx q[2];
rz(-1.3907593) q[2];
sx q[2];
rz(-1.5370739) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.20470141) q[1];
sx q[1];
rz(-1.1549779) q[1];
sx q[1];
rz(2.4067626) q[1];
rz(1.2792688) q[3];
sx q[3];
rz(-1.8390553) q[3];
sx q[3];
rz(-2.832178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5532316) q[2];
sx q[2];
rz(-1.0855805) q[2];
sx q[2];
rz(2.9300743) q[2];
rz(-1.5254947) q[3];
sx q[3];
rz(-2.1414521) q[3];
sx q[3];
rz(-2.8741527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78032988) q[0];
sx q[0];
rz(-1.4365256) q[0];
sx q[0];
rz(-2.8509675) q[0];
rz(-0.81193874) q[1];
sx q[1];
rz(-2.5135136) q[1];
sx q[1];
rz(1.5807349) q[1];
rz(-2.7079034) q[2];
sx q[2];
rz(-1.5885175) q[2];
sx q[2];
rz(1.3383404) q[2];
rz(-2.2253401) q[3];
sx q[3];
rz(-2.2753245) q[3];
sx q[3];
rz(0.38451274) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
