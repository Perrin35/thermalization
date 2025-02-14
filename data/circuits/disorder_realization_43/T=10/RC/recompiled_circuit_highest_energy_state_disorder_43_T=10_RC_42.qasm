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
rz(-2.4515406) q[0];
sx q[0];
rz(-2.2909988) q[0];
sx q[0];
rz(2.9414862) q[0];
rz(-0.19835681) q[1];
sx q[1];
rz(-0.52302066) q[1];
sx q[1];
rz(1.3318292) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057077335) q[0];
sx q[0];
rz(-0.50327089) q[0];
sx q[0];
rz(2.0664735) q[0];
rz(-pi) q[1];
rz(-0.66352377) q[2];
sx q[2];
rz(-0.71394701) q[2];
sx q[2];
rz(-2.6235142) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6864904) q[1];
sx q[1];
rz(-2.2481866) q[1];
sx q[1];
rz(-1.6325006) q[1];
rz(2.0468007) q[3];
sx q[3];
rz(-1.4972357) q[3];
sx q[3];
rz(-1.6122163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0693822) q[2];
sx q[2];
rz(-0.36082265) q[2];
sx q[2];
rz(1.8270095) q[2];
rz(-2.2383111) q[3];
sx q[3];
rz(-1.3638146) q[3];
sx q[3];
rz(-2.7956853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2717993) q[0];
sx q[0];
rz(-1.6558187) q[0];
sx q[0];
rz(-0.85579175) q[0];
rz(0.77404147) q[1];
sx q[1];
rz(-1.1639405) q[1];
sx q[1];
rz(-0.78786293) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1163569) q[0];
sx q[0];
rz(-1.5486002) q[0];
sx q[0];
rz(-1.8216013) q[0];
rz(-pi) q[1];
rz(-0.34423265) q[2];
sx q[2];
rz(-1.6107127) q[2];
sx q[2];
rz(1.0913864) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1545002) q[1];
sx q[1];
rz(-0.56114158) q[1];
sx q[1];
rz(2.2461056) q[1];
rz(-1.7807021) q[3];
sx q[3];
rz(-2.2271101) q[3];
sx q[3];
rz(2.8950952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.065980109) q[2];
sx q[2];
rz(-1.9115261) q[2];
sx q[2];
rz(2.7454929) q[2];
rz(1.5710477) q[3];
sx q[3];
rz(-1.6377662) q[3];
sx q[3];
rz(1.7861231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1876672) q[0];
sx q[0];
rz(-1.4952156) q[0];
sx q[0];
rz(-3.0809825) q[0];
rz(0.68471471) q[1];
sx q[1];
rz(-1.4002607) q[1];
sx q[1];
rz(2.8025467) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7386757) q[0];
sx q[0];
rz(-2.5581723) q[0];
sx q[0];
rz(3.1039799) q[0];
rz(2.9525063) q[2];
sx q[2];
rz(-1.4497309) q[2];
sx q[2];
rz(-2.6775286) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0319034) q[1];
sx q[1];
rz(-1.2983783) q[1];
sx q[1];
rz(1.8930356) q[1];
rz(3.1296785) q[3];
sx q[3];
rz(-0.77687009) q[3];
sx q[3];
rz(0.35038951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.48245779) q[2];
sx q[2];
rz(-1.8637916) q[2];
sx q[2];
rz(-2.6711312) q[2];
rz(-0.86236924) q[3];
sx q[3];
rz(-2.6926398) q[3];
sx q[3];
rz(-1.5615777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3622793) q[0];
sx q[0];
rz(-2.0231495) q[0];
sx q[0];
rz(0.40801868) q[0];
rz(-1.7904003) q[1];
sx q[1];
rz(-1.2974757) q[1];
sx q[1];
rz(-1.8203576) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9269598) q[0];
sx q[0];
rz(-2.0713191) q[0];
sx q[0];
rz(1.2248125) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6826019) q[2];
sx q[2];
rz(-2.4833792) q[2];
sx q[2];
rz(1.7786225) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3776748) q[1];
sx q[1];
rz(-1.6435247) q[1];
sx q[1];
rz(-3.0200028) q[1];
x q[2];
rz(0.77797555) q[3];
sx q[3];
rz(-0.56853349) q[3];
sx q[3];
rz(-1.83873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9185751) q[2];
sx q[2];
rz(-1.7097946) q[2];
sx q[2];
rz(-0.72285405) q[2];
rz(-0.84960788) q[3];
sx q[3];
rz(-1.1725715) q[3];
sx q[3];
rz(-2.1037219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1866813) q[0];
sx q[0];
rz(-0.51690042) q[0];
sx q[0];
rz(-2.5816259) q[0];
rz(2.9365183) q[1];
sx q[1];
rz(-1.2566902) q[1];
sx q[1];
rz(-2.8505039) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4093035) q[0];
sx q[0];
rz(-1.8390053) q[0];
sx q[0];
rz(1.3951673) q[0];
x q[1];
rz(0.061441378) q[2];
sx q[2];
rz(-0.80020088) q[2];
sx q[2];
rz(-2.7903008) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6711802) q[1];
sx q[1];
rz(-2.414783) q[1];
sx q[1];
rz(1.2695168) q[1];
x q[2];
rz(0.17125968) q[3];
sx q[3];
rz(-2.3823839) q[3];
sx q[3];
rz(-1.9943135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2228955) q[2];
sx q[2];
rz(-1.4783858) q[2];
sx q[2];
rz(0.29829868) q[2];
rz(-1.1693303) q[3];
sx q[3];
rz(-2.1746641) q[3];
sx q[3];
rz(-1.5112618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(2.8339612) q[0];
sx q[0];
rz(-1.6809604) q[0];
sx q[0];
rz(-0.6024012) q[0];
rz(2.7700453) q[1];
sx q[1];
rz(-1.5967775) q[1];
sx q[1];
rz(2.1098302) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96249798) q[0];
sx q[0];
rz(-2.966605) q[0];
sx q[0];
rz(1.4797158) q[0];
x q[1];
rz(1.363109) q[2];
sx q[2];
rz(-1.9676932) q[2];
sx q[2];
rz(-2.7092421) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.18373016) q[1];
sx q[1];
rz(-2.0477844) q[1];
sx q[1];
rz(1.7836003) q[1];
rz(-2.5661181) q[3];
sx q[3];
rz(-2.1361975) q[3];
sx q[3];
rz(2.7421212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3370257) q[2];
sx q[2];
rz(-2.0402543) q[2];
sx q[2];
rz(-1.5186914) q[2];
rz(2.2931781) q[3];
sx q[3];
rz(-0.87526667) q[3];
sx q[3];
rz(-3.1019822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0254211) q[0];
sx q[0];
rz(-3.0905368) q[0];
sx q[0];
rz(2.1436932) q[0];
rz(-2.8669224) q[1];
sx q[1];
rz(-0.88013595) q[1];
sx q[1];
rz(1.7707228) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8664198) q[0];
sx q[0];
rz(-0.93979561) q[0];
sx q[0];
rz(-1.5541398) q[0];
x q[1];
rz(1.9657126) q[2];
sx q[2];
rz(-0.38287258) q[2];
sx q[2];
rz(1.9667786) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.22860195) q[1];
sx q[1];
rz(-0.96768846) q[1];
sx q[1];
rz(1.8643537) q[1];
x q[2];
rz(2.3207449) q[3];
sx q[3];
rz(-2.2162262) q[3];
sx q[3];
rz(-2.5643333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2844598) q[2];
sx q[2];
rz(-0.94228116) q[2];
sx q[2];
rz(3.1119463) q[2];
rz(-2.5263785) q[3];
sx q[3];
rz(-1.7447724) q[3];
sx q[3];
rz(0.342338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-3.0520332) q[0];
sx q[0];
rz(-0.46228662) q[0];
sx q[0];
rz(-0.7830559) q[0];
rz(1.1198593) q[1];
sx q[1];
rz(-0.66711396) q[1];
sx q[1];
rz(-1.7436183) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92358855) q[0];
sx q[0];
rz(-0.62380416) q[0];
sx q[0];
rz(-2.3179834) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9080284) q[2];
sx q[2];
rz(-1.8478881) q[2];
sx q[2];
rz(-1.8943) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74136342) q[1];
sx q[1];
rz(-1.5619763) q[1];
sx q[1];
rz(-2.2386902) q[1];
rz(-0.21381883) q[3];
sx q[3];
rz(-0.79584661) q[3];
sx q[3];
rz(2.4406432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.29584259) q[2];
sx q[2];
rz(-2.5551899) q[2];
sx q[2];
rz(-1.8208549) q[2];
rz(2.1218421) q[3];
sx q[3];
rz(-1.7176065) q[3];
sx q[3];
rz(-1.7216871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9418697) q[0];
sx q[0];
rz(-2.0482735) q[0];
sx q[0];
rz(0.95440188) q[0];
rz(1.154254) q[1];
sx q[1];
rz(-2.4955165) q[1];
sx q[1];
rz(2.0571713) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1444839) q[0];
sx q[0];
rz(-2.2786744) q[0];
sx q[0];
rz(3.0860712) q[0];
x q[1];
rz(-1.2899288) q[2];
sx q[2];
rz(-1.3251588) q[2];
sx q[2];
rz(-2.7410067) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3855648) q[1];
sx q[1];
rz(-1.6431019) q[1];
sx q[1];
rz(2.0438309) q[1];
rz(0.74862759) q[3];
sx q[3];
rz(-1.4350865) q[3];
sx q[3];
rz(2.3627797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0584917) q[2];
sx q[2];
rz(-2.9497171) q[2];
sx q[2];
rz(-2.7044738) q[2];
rz(1.7982177) q[3];
sx q[3];
rz(-1.3513887) q[3];
sx q[3];
rz(3.0962211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042353543) q[0];
sx q[0];
rz(-2.1411242) q[0];
sx q[0];
rz(1.4105256) q[0];
rz(-0.22077416) q[1];
sx q[1];
rz(-0.53123728) q[1];
sx q[1];
rz(-0.9332307) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8659023) q[0];
sx q[0];
rz(-2.3161016) q[0];
sx q[0];
rz(1.5893713) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61874091) q[2];
sx q[2];
rz(-1.8685942) q[2];
sx q[2];
rz(-2.1190475) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1677959) q[1];
sx q[1];
rz(-0.78404155) q[1];
sx q[1];
rz(1.382949) q[1];
rz(-pi) q[2];
rz(-0.93818061) q[3];
sx q[3];
rz(-1.6611735) q[3];
sx q[3];
rz(0.5394906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.431388) q[2];
sx q[2];
rz(-0.34222558) q[2];
sx q[2];
rz(-0.63146511) q[2];
rz(-0.7190052) q[3];
sx q[3];
rz(-1.3471666) q[3];
sx q[3];
rz(-2.9001111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63856335) q[0];
sx q[0];
rz(-1.4894435) q[0];
sx q[0];
rz(1.3110934) q[0];
rz(2.6999264) q[1];
sx q[1];
rz(-1.3413981) q[1];
sx q[1];
rz(-1.9955019) q[1];
rz(0.031671192) q[2];
sx q[2];
rz(-1.9513672) q[2];
sx q[2];
rz(0.58917108) q[2];
rz(0.56959116) q[3];
sx q[3];
rz(-1.4841595) q[3];
sx q[3];
rz(-1.0356173) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
