OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.54685932) q[0];
sx q[0];
rz(-1.62513) q[0];
sx q[0];
rz(-0.2642785) q[0];
rz(-0.9737941) q[1];
sx q[1];
rz(5.073054) q[1];
sx q[1];
rz(10.160025) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41266325) q[0];
sx q[0];
rz(-0.67617765) q[0];
sx q[0];
rz(2.9039608) q[0];
x q[1];
rz(-2.5392883) q[2];
sx q[2];
rz(-2.3659083) q[2];
sx q[2];
rz(1.3210981) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.6890251) q[1];
sx q[1];
rz(-1.2847632) q[1];
sx q[1];
rz(0.16098117) q[1];
rz(-pi) q[2];
rz(3.0987708) q[3];
sx q[3];
rz(-0.58086568) q[3];
sx q[3];
rz(-2.7659741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.66951093) q[2];
sx q[2];
rz(-1.3005723) q[2];
sx q[2];
rz(2.0377339) q[2];
rz(-1.2708698) q[3];
sx q[3];
rz(-1.2277675) q[3];
sx q[3];
rz(0.27403533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1141777) q[0];
sx q[0];
rz(-0.52863055) q[0];
sx q[0];
rz(-2.7052178) q[0];
rz(-2.6787058) q[1];
sx q[1];
rz(-1.0375689) q[1];
sx q[1];
rz(2.8754821) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55914315) q[0];
sx q[0];
rz(-0.83280116) q[0];
sx q[0];
rz(1.0484496) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9782148) q[2];
sx q[2];
rz(-2.1401569) q[2];
sx q[2];
rz(-0.53158224) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4437342) q[1];
sx q[1];
rz(-2.2599972) q[1];
sx q[1];
rz(-0.84390784) q[1];
rz(-1.3105884) q[3];
sx q[3];
rz(-2.5187413) q[3];
sx q[3];
rz(-0.42850307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46488547) q[2];
sx q[2];
rz(-1.2767982) q[2];
sx q[2];
rz(2.6300988) q[2];
rz(-2.3320847) q[3];
sx q[3];
rz(-1.5313238) q[3];
sx q[3];
rz(-0.28765837) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4354316) q[0];
sx q[0];
rz(-1.5478739) q[0];
sx q[0];
rz(-0.92873746) q[0];
rz(-1.4061032) q[1];
sx q[1];
rz(-0.7000674) q[1];
sx q[1];
rz(1.3471289) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72235332) q[0];
sx q[0];
rz(-1.1142715) q[0];
sx q[0];
rz(-1.8629575) q[0];
x q[1];
rz(1.8560266) q[2];
sx q[2];
rz(-1.7226379) q[2];
sx q[2];
rz(1.0361995) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.24088629) q[1];
sx q[1];
rz(-0.97929231) q[1];
sx q[1];
rz(-2.790931) q[1];
rz(-pi) q[2];
rz(-1.7473162) q[3];
sx q[3];
rz(-2.4268097) q[3];
sx q[3];
rz(0.96283462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1469664) q[2];
sx q[2];
rz(-1.7549843) q[2];
sx q[2];
rz(1.6195126) q[2];
rz(0.26432031) q[3];
sx q[3];
rz(-1.0364573) q[3];
sx q[3];
rz(0.49595293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3494444) q[0];
sx q[0];
rz(-1.9932207) q[0];
sx q[0];
rz(0.96570063) q[0];
rz(2.4194338) q[1];
sx q[1];
rz(-1.637371) q[1];
sx q[1];
rz(-0.55975634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0069665466) q[0];
sx q[0];
rz(-2.7077733) q[0];
sx q[0];
rz(0.19402786) q[0];
rz(-pi) q[1];
rz(-0.83300029) q[2];
sx q[2];
rz(-0.90655316) q[2];
sx q[2];
rz(0.26495648) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8371115) q[1];
sx q[1];
rz(-1.5251535) q[1];
sx q[1];
rz(0.95506217) q[1];
rz(-2.3523931) q[3];
sx q[3];
rz(-1.7390828) q[3];
sx q[3];
rz(-2.5326953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.42797783) q[2];
sx q[2];
rz(-1.5087936) q[2];
sx q[2];
rz(-1.4245865) q[2];
rz(2.8811841) q[3];
sx q[3];
rz(-1.7491165) q[3];
sx q[3];
rz(-2.7105455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99073064) q[0];
sx q[0];
rz(-1.6435511) q[0];
sx q[0];
rz(-2.7752303) q[0];
rz(-1.5953966) q[1];
sx q[1];
rz(-2.5876744) q[1];
sx q[1];
rz(-2.7979134) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8923963) q[0];
sx q[0];
rz(-2.2049892) q[0];
sx q[0];
rz(-2.2618494) q[0];
rz(-pi) q[1];
rz(2.5630066) q[2];
sx q[2];
rz(-0.95539504) q[2];
sx q[2];
rz(-1.7631284) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1829454) q[1];
sx q[1];
rz(-1.4970386) q[1];
sx q[1];
rz(1.0720836) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.120317) q[3];
sx q[3];
rz(-2.2145503) q[3];
sx q[3];
rz(1.4847886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3917824) q[2];
sx q[2];
rz(-2.2822773) q[2];
sx q[2];
rz(3.0878477) q[2];
rz(-1.404445) q[3];
sx q[3];
rz(-2.6440547) q[3];
sx q[3];
rz(0.29156175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6546201) q[0];
sx q[0];
rz(-2.4916861) q[0];
sx q[0];
rz(-0.75575954) q[0];
rz(-3.1164363) q[1];
sx q[1];
rz(-0.92725602) q[1];
sx q[1];
rz(-2.8818534) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.320967) q[0];
sx q[0];
rz(-2.6532986) q[0];
sx q[0];
rz(0.85296209) q[0];
x q[1];
rz(1.2345215) q[2];
sx q[2];
rz(-1.5521126) q[2];
sx q[2];
rz(-0.74355723) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0745084) q[1];
sx q[1];
rz(-1.2423406) q[1];
sx q[1];
rz(-2.902608) q[1];
rz(2.0201683) q[3];
sx q[3];
rz(-1.5869889) q[3];
sx q[3];
rz(1.3390954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6529237) q[2];
sx q[2];
rz(-2.3355464) q[2];
sx q[2];
rz(3.0409813) q[2];
rz(-0.18209022) q[3];
sx q[3];
rz(-0.87825769) q[3];
sx q[3];
rz(-1.7939059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1473734) q[0];
sx q[0];
rz(-1.7213151) q[0];
sx q[0];
rz(0.31016645) q[0];
rz(2.639333) q[1];
sx q[1];
rz(-0.52400932) q[1];
sx q[1];
rz(0.60595864) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2180442) q[0];
sx q[0];
rz(-2.179562) q[0];
sx q[0];
rz(1.2952842) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.546054) q[2];
sx q[2];
rz(-2.7936802) q[2];
sx q[2];
rz(0.66875848) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.088674) q[1];
sx q[1];
rz(-2.5175736) q[1];
sx q[1];
rz(-2.1799929) q[1];
x q[2];
rz(1.0244272) q[3];
sx q[3];
rz(-1.1064648) q[3];
sx q[3];
rz(1.446561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9014152) q[2];
sx q[2];
rz(-1.1837974) q[2];
sx q[2];
rz(2.288738) q[2];
rz(-1.7715706) q[3];
sx q[3];
rz(-1.6829237) q[3];
sx q[3];
rz(0.20496932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2296427) q[0];
sx q[0];
rz(-0.59597534) q[0];
sx q[0];
rz(1.6802616) q[0];
rz(1.4029067) q[1];
sx q[1];
rz(-0.97424126) q[1];
sx q[1];
rz(3.0775552) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2647576) q[0];
sx q[0];
rz(-1.7581853) q[0];
sx q[0];
rz(0.95215709) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0766826) q[2];
sx q[2];
rz(-1.8407028) q[2];
sx q[2];
rz(0.050886521) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.774051) q[1];
sx q[1];
rz(-1.5343752) q[1];
sx q[1];
rz(-0.35518412) q[1];
rz(-3.1021032) q[3];
sx q[3];
rz(-1.2590623) q[3];
sx q[3];
rz(-1.5962275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0503851) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(-1.5554265) q[2];
rz(-2.2533916) q[3];
sx q[3];
rz(-1.9017838) q[3];
sx q[3];
rz(-0.92938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14426194) q[0];
sx q[0];
rz(-2.0781131) q[0];
sx q[0];
rz(-1.4319179) q[0];
rz(0.56888467) q[1];
sx q[1];
rz(-2.6061997) q[1];
sx q[1];
rz(-2.0137537) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9147946) q[0];
sx q[0];
rz(-2.9920122) q[0];
sx q[0];
rz(0.10766715) q[0];
x q[1];
rz(2.2531829) q[2];
sx q[2];
rz(-2.3896304) q[2];
sx q[2];
rz(1.0351406) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9319218) q[1];
sx q[1];
rz(-1.5591963) q[1];
sx q[1];
rz(1.8124707) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3663775) q[3];
sx q[3];
rz(-1.0421703) q[3];
sx q[3];
rz(-0.60929326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.613712) q[2];
sx q[2];
rz(-0.85313672) q[2];
sx q[2];
rz(-0.17364994) q[2];
rz(0.33637834) q[3];
sx q[3];
rz(-1.222638) q[3];
sx q[3];
rz(-0.39150795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7062475) q[0];
sx q[0];
rz(-2.5543537) q[0];
sx q[0];
rz(-1.6760814) q[0];
rz(-2.3174875) q[1];
sx q[1];
rz(-1.6128287) q[1];
sx q[1];
rz(-2.5691659) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7596282) q[0];
sx q[0];
rz(-1.0257162) q[0];
sx q[0];
rz(2.1419924) q[0];
rz(-0.27406759) q[2];
sx q[2];
rz(-1.1268508) q[2];
sx q[2];
rz(2.519671) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.010667) q[1];
sx q[1];
rz(-1.3548684) q[1];
sx q[1];
rz(0.24433498) q[1];
rz(-pi) q[2];
rz(1.8626067) q[3];
sx q[3];
rz(-2.2581873) q[3];
sx q[3];
rz(-0.22112267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3616025) q[2];
sx q[2];
rz(-0.47452351) q[2];
sx q[2];
rz(-0.53722107) q[2];
rz(-1.0572664) q[3];
sx q[3];
rz(-0.89151645) q[3];
sx q[3];
rz(0.69361544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2873516) q[0];
sx q[0];
rz(-1.1653405) q[0];
sx q[0];
rz(-1.5821138) q[0];
rz(-0.46335012) q[1];
sx q[1];
rz(-0.87711038) q[1];
sx q[1];
rz(-1.6323485) q[1];
rz(-2.2255185) q[2];
sx q[2];
rz(-1.6486042) q[2];
sx q[2];
rz(-0.83124607) q[2];
rz(1.0827071) q[3];
sx q[3];
rz(-1.6255717) q[3];
sx q[3];
rz(-2.0902904) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];