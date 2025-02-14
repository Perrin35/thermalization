OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.4888157) q[0];
sx q[0];
rz(-0.96537679) q[0];
sx q[0];
rz(2.2247347) q[0];
rz(1.4404453) q[1];
sx q[1];
rz(-1.0280161) q[1];
sx q[1];
rz(2.4457959) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0442565) q[0];
sx q[0];
rz(-1.6249661) q[0];
sx q[0];
rz(-1.9977117) q[0];
x q[1];
rz(-0.037362413) q[2];
sx q[2];
rz(-1.6184095) q[2];
sx q[2];
rz(-1.3726774) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2604044) q[1];
sx q[1];
rz(-1.9450448) q[1];
sx q[1];
rz(1.7541691) q[1];
x q[2];
rz(0.37000044) q[3];
sx q[3];
rz(-1.1048855) q[3];
sx q[3];
rz(0.53088684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9838788) q[2];
sx q[2];
rz(-0.82008755) q[2];
sx q[2];
rz(2.2411818) q[2];
rz(2.3484717) q[3];
sx q[3];
rz(-1.4916865) q[3];
sx q[3];
rz(0.69860631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6124436) q[0];
sx q[0];
rz(-1.9731033) q[0];
sx q[0];
rz(2.4560665) q[0];
rz(0.43288747) q[1];
sx q[1];
rz(-1.2665749) q[1];
sx q[1];
rz(1.8276851) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1788119) q[0];
sx q[0];
rz(-1.0355486) q[0];
sx q[0];
rz(-0.58387941) q[0];
rz(-pi) q[1];
rz(2.6816112) q[2];
sx q[2];
rz(-2.26074) q[2];
sx q[2];
rz(-1.3342109) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.59456277) q[1];
sx q[1];
rz(-2.4395151) q[1];
sx q[1];
rz(-2.3365993) q[1];
x q[2];
rz(0.39137381) q[3];
sx q[3];
rz(-0.59283601) q[3];
sx q[3];
rz(-0.22030574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1757397) q[2];
sx q[2];
rz(-2.7312835) q[2];
sx q[2];
rz(0.24743323) q[2];
rz(-3.0257709) q[3];
sx q[3];
rz(-1.7146866) q[3];
sx q[3];
rz(0.57870948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9847357) q[0];
sx q[0];
rz(-0.37608376) q[0];
sx q[0];
rz(1.0648741) q[0];
rz(2.5085326) q[1];
sx q[1];
rz(-1.1587016) q[1];
sx q[1];
rz(3.0050468) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8234607) q[0];
sx q[0];
rz(-0.7366418) q[0];
sx q[0];
rz(0.036286906) q[0];
x q[1];
rz(-1.1077131) q[2];
sx q[2];
rz(-1.255724) q[2];
sx q[2];
rz(1.7618432) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8571843) q[1];
sx q[1];
rz(-2.2940002) q[1];
sx q[1];
rz(-0.64874362) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7259215) q[3];
sx q[3];
rz(-1.1729673) q[3];
sx q[3];
rz(-1.3597387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.69918767) q[2];
sx q[2];
rz(-1.8892611) q[2];
sx q[2];
rz(1.5415972) q[2];
rz(2.048118) q[3];
sx q[3];
rz(-1.6396061) q[3];
sx q[3];
rz(-2.7961965) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4271451) q[0];
sx q[0];
rz(-2.8801425) q[0];
sx q[0];
rz(-0.57681042) q[0];
rz(-2.414074) q[1];
sx q[1];
rz(-1.0914165) q[1];
sx q[1];
rz(2.2732546) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.627358) q[0];
sx q[0];
rz(-1.9666202) q[0];
sx q[0];
rz(2.4577599) q[0];
rz(2.6170591) q[2];
sx q[2];
rz(-1.0528334) q[2];
sx q[2];
rz(0.7690767) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7817745) q[1];
sx q[1];
rz(-1.9137772) q[1];
sx q[1];
rz(0.69056679) q[1];
x q[2];
rz(-2.6336977) q[3];
sx q[3];
rz(-2.1334071) q[3];
sx q[3];
rz(-0.25319448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9463828) q[2];
sx q[2];
rz(-1.3888487) q[2];
sx q[2];
rz(-0.69551224) q[2];
rz(-0.21466151) q[3];
sx q[3];
rz(-1.6404459) q[3];
sx q[3];
rz(-2.3282839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6014086) q[0];
sx q[0];
rz(-2.5267127) q[0];
sx q[0];
rz(1.1001128) q[0];
rz(-0.21057883) q[1];
sx q[1];
rz(-0.63340488) q[1];
sx q[1];
rz(1.6483773) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57015145) q[0];
sx q[0];
rz(-2.049027) q[0];
sx q[0];
rz(-0.5334665) q[0];
rz(-pi) q[1];
rz(1.3927685) q[2];
sx q[2];
rz(-1.4881987) q[2];
sx q[2];
rz(2.1796153) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4752057) q[1];
sx q[1];
rz(-2.0568486) q[1];
sx q[1];
rz(0.17611019) q[1];
rz(-pi) q[2];
rz(1.8214956) q[3];
sx q[3];
rz(-2.2226102) q[3];
sx q[3];
rz(-2.1194315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7275927) q[2];
sx q[2];
rz(-1.4294581) q[2];
sx q[2];
rz(0.89699888) q[2];
rz(-1.2796848) q[3];
sx q[3];
rz(-1.0106267) q[3];
sx q[3];
rz(0.047209386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9996662) q[0];
sx q[0];
rz(-2.5157295) q[0];
sx q[0];
rz(2.8114787) q[0];
rz(-0.64078981) q[1];
sx q[1];
rz(-1.7658486) q[1];
sx q[1];
rz(0.97553387) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3307083) q[0];
sx q[0];
rz(-2.160794) q[0];
sx q[0];
rz(1.6390653) q[0];
x q[1];
rz(0.26872608) q[2];
sx q[2];
rz(-1.8789504) q[2];
sx q[2];
rz(1.9374963) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9258049) q[1];
sx q[1];
rz(-2.6138209) q[1];
sx q[1];
rz(2.5592519) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1485432) q[3];
sx q[3];
rz(-1.7429757) q[3];
sx q[3];
rz(0.14587054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2834125) q[2];
sx q[2];
rz(-2.8874669) q[2];
sx q[2];
rz(0.88968366) q[2];
rz(0.59600082) q[3];
sx q[3];
rz(-2.0723876) q[3];
sx q[3];
rz(3.0349351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0058873) q[0];
sx q[0];
rz(-1.1361253) q[0];
sx q[0];
rz(2.0086052) q[0];
rz(0.63944447) q[1];
sx q[1];
rz(-2.0638128) q[1];
sx q[1];
rz(2.1741672) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15617642) q[0];
sx q[0];
rz(-1.7207017) q[0];
sx q[0];
rz(2.085859) q[0];
rz(2.7616639) q[2];
sx q[2];
rz(-1.8779781) q[2];
sx q[2];
rz(1.2387222) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.51951197) q[1];
sx q[1];
rz(-1.4556307) q[1];
sx q[1];
rz(0.22611109) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2626454) q[3];
sx q[3];
rz(-2.2053092) q[3];
sx q[3];
rz(0.18116118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7449259) q[2];
sx q[2];
rz(-2.1325839) q[2];
sx q[2];
rz(0.94433707) q[2];
rz(2.5189597) q[3];
sx q[3];
rz(-2.1532652) q[3];
sx q[3];
rz(-2.3545789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(3.0866163) q[0];
sx q[0];
rz(-0.80626196) q[0];
sx q[0];
rz(-0.0042313519) q[0];
rz(-1.6731693) q[1];
sx q[1];
rz(-2.3783042) q[1];
sx q[1];
rz(0.17446336) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8426822) q[0];
sx q[0];
rz(-1.7571661) q[0];
sx q[0];
rz(3.0475265) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44850839) q[2];
sx q[2];
rz(-1.4786198) q[2];
sx q[2];
rz(-1.1223457) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9618865) q[1];
sx q[1];
rz(-0.70155662) q[1];
sx q[1];
rz(1.6732193) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8352083) q[3];
sx q[3];
rz(-0.2443265) q[3];
sx q[3];
rz(0.78247386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8909495) q[2];
sx q[2];
rz(-3.0477016) q[2];
sx q[2];
rz(-0.54035652) q[2];
rz(-3.12449) q[3];
sx q[3];
rz(-0.98186866) q[3];
sx q[3];
rz(0.66177773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0509725) q[0];
sx q[0];
rz(-2.3391188) q[0];
sx q[0];
rz(-0.78659868) q[0];
rz(-2.0358918) q[1];
sx q[1];
rz(-1.4812171) q[1];
sx q[1];
rz(0.22122637) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1329256) q[0];
sx q[0];
rz(-1.8485539) q[0];
sx q[0];
rz(-0.8322258) q[0];
x q[1];
rz(1.073313) q[2];
sx q[2];
rz(-0.19693298) q[2];
sx q[2];
rz(-0.54468583) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7936642) q[1];
sx q[1];
rz(-1.3634487) q[1];
sx q[1];
rz(0.4936895) q[1];
rz(1.33696) q[3];
sx q[3];
rz(-0.89409308) q[3];
sx q[3];
rz(1.9399583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1456387) q[2];
sx q[2];
rz(-0.24253878) q[2];
sx q[2];
rz(1.185574) q[2];
rz(1.1061741) q[3];
sx q[3];
rz(-2.1198876) q[3];
sx q[3];
rz(-0.98418981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36158654) q[0];
sx q[0];
rz(-1.3573283) q[0];
sx q[0];
rz(-0.7097882) q[0];
rz(0.89372006) q[1];
sx q[1];
rz(-0.62785405) q[1];
sx q[1];
rz(-0.46700221) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92999536) q[0];
sx q[0];
rz(-2.8401137) q[0];
sx q[0];
rz(1.0903574) q[0];
rz(-pi) q[1];
rz(0.94552083) q[2];
sx q[2];
rz(-1.7904803) q[2];
sx q[2];
rz(-1.4131119) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.71605395) q[1];
sx q[1];
rz(-1.8093093) q[1];
sx q[1];
rz(3.1032326) q[1];
rz(-0.8740467) q[3];
sx q[3];
rz(-2.120904) q[3];
sx q[3];
rz(2.3693905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3422157) q[2];
sx q[2];
rz(-0.90505427) q[2];
sx q[2];
rz(-3.0055935) q[2];
rz(2.811725) q[3];
sx q[3];
rz(-2.0743275) q[3];
sx q[3];
rz(-0.29712591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8316523) q[0];
sx q[0];
rz(-1.1158054) q[0];
sx q[0];
rz(0.62624395) q[0];
rz(-0.88824875) q[1];
sx q[1];
rz(-1.8969957) q[1];
sx q[1];
rz(-3.0096164) q[1];
rz(1.7952948) q[2];
sx q[2];
rz(-1.9585051) q[2];
sx q[2];
rz(-0.060123882) q[2];
rz(-2.6693564) q[3];
sx q[3];
rz(-2.355081) q[3];
sx q[3];
rz(-0.045674617) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
