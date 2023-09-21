OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3044843) q[0];
sx q[0];
rz(-1.6882856) q[0];
sx q[0];
rz(2.8300571) q[0];
rz(2.7040634) q[1];
sx q[1];
rz(-1.3181926) q[1];
sx q[1];
rz(-0.55895609) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9851202) q[0];
sx q[0];
rz(-0.44056842) q[0];
sx q[0];
rz(1.23929) q[0];
rz(-pi) q[1];
rz(-0.55732255) q[2];
sx q[2];
rz(-1.5601336) q[2];
sx q[2];
rz(-2.9023841) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6701339) q[1];
sx q[1];
rz(-2.11073) q[1];
sx q[1];
rz(2.259841) q[1];
rz(-2.7028014) q[3];
sx q[3];
rz(-1.8306797) q[3];
sx q[3];
rz(-2.4364542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3922334) q[2];
sx q[2];
rz(-1.2831251) q[2];
sx q[2];
rz(-2.5048845) q[2];
rz(-2.2926245) q[3];
sx q[3];
rz(-2.520112) q[3];
sx q[3];
rz(-0.23392114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37671509) q[0];
sx q[0];
rz(-2.8945518) q[0];
sx q[0];
rz(0.15287457) q[0];
rz(-2.3846467) q[1];
sx q[1];
rz(-1.5870973) q[1];
sx q[1];
rz(2.1551932) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9372285) q[0];
sx q[0];
rz(-1.6305271) q[0];
sx q[0];
rz(1.6110957) q[0];
rz(-pi) q[1];
rz(-2.5949391) q[2];
sx q[2];
rz(-0.84178998) q[2];
sx q[2];
rz(2.5559049) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4437372) q[1];
sx q[1];
rz(-1.9471696) q[1];
sx q[1];
rz(1.3156375) q[1];
rz(0.13568474) q[3];
sx q[3];
rz(-1.8968582) q[3];
sx q[3];
rz(-0.56860926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5793005) q[2];
sx q[2];
rz(-1.219517) q[2];
sx q[2];
rz(0.78312773) q[2];
rz(-0.018571818) q[3];
sx q[3];
rz(-1.637807) q[3];
sx q[3];
rz(0.40772453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0884393) q[0];
sx q[0];
rz(-2.8494371) q[0];
sx q[0];
rz(0.96167481) q[0];
rz(2.7812474) q[1];
sx q[1];
rz(-1.1018437) q[1];
sx q[1];
rz(0.12869421) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6260687) q[0];
sx q[0];
rz(-1.4883853) q[0];
sx q[0];
rz(3.0936196) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.53993291) q[2];
sx q[2];
rz(-0.88368249) q[2];
sx q[2];
rz(2.5643258) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.80850959) q[1];
sx q[1];
rz(-1.8500449) q[1];
sx q[1];
rz(0.46344325) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4314753) q[3];
sx q[3];
rz(-2.3686667) q[3];
sx q[3];
rz(-1.3611925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0814357) q[2];
sx q[2];
rz(-1.1971985) q[2];
sx q[2];
rz(-1.8998247) q[2];
rz(0.5870108) q[3];
sx q[3];
rz(-2.185052) q[3];
sx q[3];
rz(2.164042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.509165) q[0];
sx q[0];
rz(-2.2594663) q[0];
sx q[0];
rz(2.0571016) q[0];
rz(1.658461) q[1];
sx q[1];
rz(-2.5741534) q[1];
sx q[1];
rz(3.049057) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1922069) q[0];
sx q[0];
rz(-1.8398251) q[0];
sx q[0];
rz(-1.8640679) q[0];
rz(-pi) q[1];
rz(-2.7961568) q[2];
sx q[2];
rz(-2.0256809) q[2];
sx q[2];
rz(2.5330184) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9738237) q[1];
sx q[1];
rz(-1.9134221) q[1];
sx q[1];
rz(1.6325566) q[1];
x q[2];
rz(-1.362364) q[3];
sx q[3];
rz(-0.66550335) q[3];
sx q[3];
rz(3.0968551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3080421) q[2];
sx q[2];
rz(-1.4414859) q[2];
sx q[2];
rz(2.8095424) q[2];
rz(-2.0856693) q[3];
sx q[3];
rz(-2.8639586) q[3];
sx q[3];
rz(-2.5312996) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8191391) q[0];
sx q[0];
rz(-1.2912913) q[0];
sx q[0];
rz(-2.9300368) q[0];
rz(-1.8353204) q[1];
sx q[1];
rz(-1.2439367) q[1];
sx q[1];
rz(-0.64770118) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0603795) q[0];
sx q[0];
rz(-0.19506422) q[0];
sx q[0];
rz(2.5293406) q[0];
rz(1.7251882) q[2];
sx q[2];
rz(-1.9152181) q[2];
sx q[2];
rz(1.7248578) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9808637) q[1];
sx q[1];
rz(-2.594922) q[1];
sx q[1];
rz(-1.9561808) q[1];
rz(-0.6661617) q[3];
sx q[3];
rz(-0.024346711) q[3];
sx q[3];
rz(1.419988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5806879) q[2];
sx q[2];
rz(-0.40955341) q[2];
sx q[2];
rz(0.69331759) q[2];
rz(-0.66926113) q[3];
sx q[3];
rz(-1.7047434) q[3];
sx q[3];
rz(0.0049237331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82419056) q[0];
sx q[0];
rz(-0.22739246) q[0];
sx q[0];
rz(1.90907) q[0];
rz(2.0690074) q[1];
sx q[1];
rz(-2.0697846) q[1];
sx q[1];
rz(0.17428621) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.106819) q[0];
sx q[0];
rz(-2.3475721) q[0];
sx q[0];
rz(1.130571) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86485483) q[2];
sx q[2];
rz(-2.1055429) q[2];
sx q[2];
rz(0.45644444) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.41387687) q[1];
sx q[1];
rz(-1.1912279) q[1];
sx q[1];
rz(0.80645251) q[1];
rz(-pi) q[2];
rz(-1.2586081) q[3];
sx q[3];
rz(-1.6815261) q[3];
sx q[3];
rz(0.98715106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8217414) q[2];
sx q[2];
rz(-2.3345626) q[2];
sx q[2];
rz(-2.9439587) q[2];
rz(-0.28891426) q[3];
sx q[3];
rz(-2.2646326) q[3];
sx q[3];
rz(-1.4060085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.063868) q[0];
sx q[0];
rz(-0.48982319) q[0];
sx q[0];
rz(-2.9329964) q[0];
rz(0.96616191) q[1];
sx q[1];
rz(-1.1599133) q[1];
sx q[1];
rz(-1.6360412) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2710072) q[0];
sx q[0];
rz(-2.2858372) q[0];
sx q[0];
rz(-0.72512759) q[0];
x q[1];
rz(0.079300785) q[2];
sx q[2];
rz(-2.7100025) q[2];
sx q[2];
rz(1.1232131) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.70430763) q[1];
sx q[1];
rz(-2.0211126) q[1];
sx q[1];
rz(-0.98547658) q[1];
x q[2];
rz(-2.8553477) q[3];
sx q[3];
rz(-1.4713333) q[3];
sx q[3];
rz(-0.092982987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2839526) q[2];
sx q[2];
rz(-2.3985034) q[2];
sx q[2];
rz(-3.0440142) q[2];
rz(1.3939259) q[3];
sx q[3];
rz(-1.3503617) q[3];
sx q[3];
rz(0.71684366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7512648) q[0];
sx q[0];
rz(-1.8126235) q[0];
sx q[0];
rz(0.51399291) q[0];
rz(0.12318525) q[1];
sx q[1];
rz(-2.8942278) q[1];
sx q[1];
rz(-0.93200144) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4829464) q[0];
sx q[0];
rz(-2.6771149) q[0];
sx q[0];
rz(1.8770201) q[0];
rz(-2.2064662) q[2];
sx q[2];
rz(-1.7889651) q[2];
sx q[2];
rz(-0.93014923) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2145558) q[1];
sx q[1];
rz(-0.80074691) q[1];
sx q[1];
rz(2.3826249) q[1];
x q[2];
rz(-1.7979513) q[3];
sx q[3];
rz(-2.4121768) q[3];
sx q[3];
rz(-1.9279355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1121858) q[2];
sx q[2];
rz(-1.1108578) q[2];
sx q[2];
rz(0.4294447) q[2];
rz(1.2094234) q[3];
sx q[3];
rz(-0.37241396) q[3];
sx q[3];
rz(-0.69303524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3747303) q[0];
sx q[0];
rz(-1.466789) q[0];
sx q[0];
rz(1.8027579) q[0];
rz(-0.70612899) q[1];
sx q[1];
rz(-1.8920205) q[1];
sx q[1];
rz(2.1910117) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3830519) q[0];
sx q[0];
rz(-2.0630815) q[0];
sx q[0];
rz(2.5581215) q[0];
rz(-pi) q[1];
rz(1.1549994) q[2];
sx q[2];
rz(-1.6396513) q[2];
sx q[2];
rz(-1.6378251) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7710167) q[1];
sx q[1];
rz(-2.5713213) q[1];
sx q[1];
rz(-1.9221406) q[1];
rz(-pi) q[2];
rz(-1.9825963) q[3];
sx q[3];
rz(-1.3367062) q[3];
sx q[3];
rz(1.4300508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6616228) q[2];
sx q[2];
rz(-1.491549) q[2];
sx q[2];
rz(2.6573112) q[2];
rz(-2.2144923) q[3];
sx q[3];
rz(-1.8323332) q[3];
sx q[3];
rz(2.5203729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1442239) q[0];
sx q[0];
rz(-0.08865083) q[0];
sx q[0];
rz(0.22928672) q[0];
rz(0.43481049) q[1];
sx q[1];
rz(-1.2229342) q[1];
sx q[1];
rz(2.4226709) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15764788) q[0];
sx q[0];
rz(-0.63417182) q[0];
sx q[0];
rz(-1.4112524) q[0];
rz(3.0232593) q[2];
sx q[2];
rz(-0.98630691) q[2];
sx q[2];
rz(-3.0362533) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.194866) q[1];
sx q[1];
rz(-0.39992878) q[1];
sx q[1];
rz(0.3868133) q[1];
rz(-pi) q[2];
rz(0.50456725) q[3];
sx q[3];
rz(-2.1310398) q[3];
sx q[3];
rz(1.954078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3832613) q[2];
sx q[2];
rz(-1.9394082) q[2];
sx q[2];
rz(0.74679217) q[2];
rz(-2.2693999) q[3];
sx q[3];
rz(-0.82834297) q[3];
sx q[3];
rz(0.94223589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.223021) q[0];
sx q[0];
rz(-1.3970319) q[0];
sx q[0];
rz(-1.3402517) q[0];
rz(0.37721286) q[1];
sx q[1];
rz(-1.6709534) q[1];
sx q[1];
rz(0.51660641) q[1];
rz(-0.51858356) q[2];
sx q[2];
rz(-1.2972144) q[2];
sx q[2];
rz(-1.5757061) q[2];
rz(0.71729284) q[3];
sx q[3];
rz(-1.4581231) q[3];
sx q[3];
rz(0.067618528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
