OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.9956545) q[0];
sx q[0];
rz(-0.50322682) q[0];
sx q[0];
rz(-0.72416645) q[0];
rz(-2.5016298) q[1];
sx q[1];
rz(-2.6115186) q[1];
sx q[1];
rz(0.78483265) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35207507) q[0];
sx q[0];
rz(-1.2791469) q[0];
sx q[0];
rz(-1.7910936) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6720812) q[2];
sx q[2];
rz(-1.9874319) q[2];
sx q[2];
rz(-3.0243304) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9468294) q[1];
sx q[1];
rz(-2.0954872) q[1];
sx q[1];
rz(2.8835433) q[1];
rz(-pi) q[2];
rz(-1.4943487) q[3];
sx q[3];
rz(-2.7259698) q[3];
sx q[3];
rz(1.7474183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5518387) q[2];
sx q[2];
rz(-1.7244312) q[2];
sx q[2];
rz(-3.0736249) q[2];
rz(3.0170278) q[3];
sx q[3];
rz(-0.3228651) q[3];
sx q[3];
rz(1.7547866) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2215866) q[0];
sx q[0];
rz(-3.0060372) q[0];
sx q[0];
rz(0.24366972) q[0];
rz(2.5098353) q[1];
sx q[1];
rz(-1.4032204) q[1];
sx q[1];
rz(1.3557281) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6512017) q[0];
sx q[0];
rz(-0.25679195) q[0];
sx q[0];
rz(-1.6780361) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7849814) q[2];
sx q[2];
rz(-1.9372254) q[2];
sx q[2];
rz(2.8474142) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1214952) q[1];
sx q[1];
rz(-2.0165682) q[1];
sx q[1];
rz(2.9989468) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0725067) q[3];
sx q[3];
rz(-2.54832) q[3];
sx q[3];
rz(3.1309576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0791066) q[2];
sx q[2];
rz(-2.1495543) q[2];
sx q[2];
rz(2.8919354) q[2];
rz(2.6349973) q[3];
sx q[3];
rz(-1.6258312) q[3];
sx q[3];
rz(2.8095968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8963985) q[0];
sx q[0];
rz(-1.9165374) q[0];
sx q[0];
rz(-0.89843345) q[0];
rz(-1.3348745) q[1];
sx q[1];
rz(-1.2355665) q[1];
sx q[1];
rz(1.2737087) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1153591) q[0];
sx q[0];
rz(-0.44239487) q[0];
sx q[0];
rz(2.1917403) q[0];
rz(-pi) q[1];
rz(-1.6689698) q[2];
sx q[2];
rz(-1.2766826) q[2];
sx q[2];
rz(1.3054747) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9540005) q[1];
sx q[1];
rz(-2.2224269) q[1];
sx q[1];
rz(-2.0152394) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5898315) q[3];
sx q[3];
rz(-0.80358395) q[3];
sx q[3];
rz(1.7337652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0818103) q[2];
sx q[2];
rz(-0.60478294) q[2];
sx q[2];
rz(0.95345062) q[2];
rz(-3.1070784) q[3];
sx q[3];
rz(-2.3551066) q[3];
sx q[3];
rz(-0.22687337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2801441) q[0];
sx q[0];
rz(-0.21629688) q[0];
sx q[0];
rz(-0.24818534) q[0];
rz(-1.0379627) q[1];
sx q[1];
rz(-1.1231517) q[1];
sx q[1];
rz(3.0674556) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.603133) q[0];
sx q[0];
rz(-0.54766253) q[0];
sx q[0];
rz(2.0752226) q[0];
rz(-2.7934302) q[2];
sx q[2];
rz(-0.92218883) q[2];
sx q[2];
rz(-1.6096889) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3809507) q[1];
sx q[1];
rz(-1.6886854) q[1];
sx q[1];
rz(-1.2537434) q[1];
rz(0.24012633) q[3];
sx q[3];
rz(-1.4014981) q[3];
sx q[3];
rz(-1.5907767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8903824) q[2];
sx q[2];
rz(-0.40428287) q[2];
sx q[2];
rz(-3.1029491) q[2];
rz(-0.97366992) q[3];
sx q[3];
rz(-0.49574167) q[3];
sx q[3];
rz(-0.27004778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53428179) q[0];
sx q[0];
rz(-1.5357635) q[0];
sx q[0];
rz(-1.3624396) q[0];
rz(-2.3249987) q[1];
sx q[1];
rz(-1.2885619) q[1];
sx q[1];
rz(-1.978925) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8602596) q[0];
sx q[0];
rz(-1.6006002) q[0];
sx q[0];
rz(-1.6822862) q[0];
rz(-pi) q[1];
rz(-1.2478998) q[2];
sx q[2];
rz(-1.4100417) q[2];
sx q[2];
rz(1.2231183) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8414383) q[1];
sx q[1];
rz(-1.9342124) q[1];
sx q[1];
rz(2.9363948) q[1];
rz(-pi) q[2];
rz(-0.88697042) q[3];
sx q[3];
rz(-2.5081222) q[3];
sx q[3];
rz(-2.5630643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6267307) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(0.14222063) q[2];
rz(-0.90406117) q[3];
sx q[3];
rz(-1.3198493) q[3];
sx q[3];
rz(-2.952125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0267923) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(-1.4676771) q[0];
rz(-2.5698075) q[1];
sx q[1];
rz(-0.3586868) q[1];
sx q[1];
rz(0.30803672) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.98826) q[0];
sx q[0];
rz(-1.6969661) q[0];
sx q[0];
rz(1.4754962) q[0];
x q[1];
rz(0.61200895) q[2];
sx q[2];
rz(-2.0888121) q[2];
sx q[2];
rz(0.81105622) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.249237) q[1];
sx q[1];
rz(-2.4294469) q[1];
sx q[1];
rz(1.874079) q[1];
rz(-pi) q[2];
rz(2.7932348) q[3];
sx q[3];
rz(-1.6568686) q[3];
sx q[3];
rz(-1.100988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3623111) q[2];
sx q[2];
rz(-1.4004536) q[2];
sx q[2];
rz(1.9936838) q[2];
rz(2.4273196) q[3];
sx q[3];
rz(-0.89759421) q[3];
sx q[3];
rz(-0.64546293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4797453) q[0];
sx q[0];
rz(-0.84091887) q[0];
sx q[0];
rz(-3.0116144) q[0];
rz(3.1107483) q[1];
sx q[1];
rz(-1.2896616) q[1];
sx q[1];
rz(0.67108363) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13043159) q[0];
sx q[0];
rz(-1.6058308) q[0];
sx q[0];
rz(3.0879211) q[0];
rz(1.7332156) q[2];
sx q[2];
rz(-0.186609) q[2];
sx q[2];
rz(2.9090372) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.259686) q[1];
sx q[1];
rz(-2.7555008) q[1];
sx q[1];
rz(-1.2950456) q[1];
rz(3.0942261) q[3];
sx q[3];
rz(-2.2455375) q[3];
sx q[3];
rz(-2.588152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0188296) q[2];
sx q[2];
rz(-1.6100223) q[2];
sx q[2];
rz(-0.57787952) q[2];
rz(-3.1130062) q[3];
sx q[3];
rz(-1.8609906) q[3];
sx q[3];
rz(1.8813429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9776483) q[0];
sx q[0];
rz(-0.90286911) q[0];
sx q[0];
rz(0.41241616) q[0];
rz(1.4498129) q[1];
sx q[1];
rz(-1.7990566) q[1];
sx q[1];
rz(-1.9746045) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1056571) q[0];
sx q[0];
rz(-1.3577537) q[0];
sx q[0];
rz(0.94353326) q[0];
rz(-pi) q[1];
rz(0.86161676) q[2];
sx q[2];
rz(-1.3590727) q[2];
sx q[2];
rz(2.784563) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3762714) q[1];
sx q[1];
rz(-1.6011392) q[1];
sx q[1];
rz(-2.9529851) q[1];
x q[2];
rz(1.0409045) q[3];
sx q[3];
rz(-1.7003254) q[3];
sx q[3];
rz(-2.0199752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2010487) q[2];
sx q[2];
rz(-2.2622435) q[2];
sx q[2];
rz(2.9525625) q[2];
rz(-2.9947301) q[3];
sx q[3];
rz(-2.9569914) q[3];
sx q[3];
rz(1.7485025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1259595) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(0.92700672) q[0];
rz(1.758763) q[1];
sx q[1];
rz(-2.5320876) q[1];
sx q[1];
rz(-1.4896726) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4091858) q[0];
sx q[0];
rz(-2.272185) q[0];
sx q[0];
rz(-2.5138445) q[0];
rz(-2.0428558) q[2];
sx q[2];
rz(-1.5257116) q[2];
sx q[2];
rz(0.91143196) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8330194) q[1];
sx q[1];
rz(-2.8816869) q[1];
sx q[1];
rz(0.72867568) q[1];
rz(-2.3143523) q[3];
sx q[3];
rz(-0.67875553) q[3];
sx q[3];
rz(0.37952207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5513409) q[2];
sx q[2];
rz(-0.44289032) q[2];
sx q[2];
rz(1.4302953) q[2];
rz(0.57724214) q[3];
sx q[3];
rz(-0.8876628) q[3];
sx q[3];
rz(1.757471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5883314) q[0];
sx q[0];
rz(-1.3681148) q[0];
sx q[0];
rz(-2.8531895) q[0];
rz(2.6092031) q[1];
sx q[1];
rz(-0.45982292) q[1];
sx q[1];
rz(2.9945701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4924016) q[0];
sx q[0];
rz(-0.99897879) q[0];
sx q[0];
rz(-0.4581106) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9016501) q[2];
sx q[2];
rz(-1.502617) q[2];
sx q[2];
rz(2.6477674) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3210996) q[1];
sx q[1];
rz(-1.5309146) q[1];
sx q[1];
rz(-2.189373) q[1];
rz(0.57659984) q[3];
sx q[3];
rz(-1.4398265) q[3];
sx q[3];
rz(0.031158202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0816575) q[2];
sx q[2];
rz(-0.43490484) q[2];
sx q[2];
rz(2.3975513) q[2];
rz(-0.75731164) q[3];
sx q[3];
rz(-1.7777187) q[3];
sx q[3];
rz(-1.7448759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.025678) q[0];
sx q[0];
rz(-2.0712576) q[0];
sx q[0];
rz(2.0448137) q[0];
rz(0.81746447) q[1];
sx q[1];
rz(-1.9348963) q[1];
sx q[1];
rz(2.5111326) q[1];
rz(0.4521162) q[2];
sx q[2];
rz(-1.631626) q[2];
sx q[2];
rz(0.22125868) q[2];
rz(1.3656473) q[3];
sx q[3];
rz(-2.57006) q[3];
sx q[3];
rz(-0.80001696) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
