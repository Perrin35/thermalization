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
rz(1.9873729) q[0];
sx q[0];
rz(4.6648751) q[0];
sx q[0];
rz(9.2821791) q[0];
rz(-2.5481186) q[1];
sx q[1];
rz(-0.65294099) q[1];
sx q[1];
rz(-2.0963734) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3680688) q[0];
sx q[0];
rz(-1.1127825) q[0];
sx q[0];
rz(2.5450443) q[0];
x q[1];
rz(-1.6204616) q[2];
sx q[2];
rz(-0.6383183) q[2];
sx q[2];
rz(2.3172003) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3821311) q[1];
sx q[1];
rz(-2.9540017) q[1];
sx q[1];
rz(-2.1407024) q[1];
x q[2];
rz(1.0020489) q[3];
sx q[3];
rz(-1.6709329) q[3];
sx q[3];
rz(1.9619693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.97592252) q[2];
sx q[2];
rz(-2.133635) q[2];
sx q[2];
rz(-2.9284076) q[2];
rz(0.91607696) q[3];
sx q[3];
rz(-0.35917425) q[3];
sx q[3];
rz(-0.38755125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85584545) q[0];
sx q[0];
rz(-2.2540932) q[0];
sx q[0];
rz(-2.6708653) q[0];
rz(1.1031021) q[1];
sx q[1];
rz(-0.50869894) q[1];
sx q[1];
rz(-2.5618166) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8743917) q[0];
sx q[0];
rz(-1.6306595) q[0];
sx q[0];
rz(0.52459985) q[0];
x q[1];
rz(-2.5986541) q[2];
sx q[2];
rz(-1.8576872) q[2];
sx q[2];
rz(2.8101588) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3892711) q[1];
sx q[1];
rz(-1.0022517) q[1];
sx q[1];
rz(-2.6804377) q[1];
rz(-1.2006496) q[3];
sx q[3];
rz(-1.8164299) q[3];
sx q[3];
rz(-2.4963238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1595999) q[2];
sx q[2];
rz(-0.85593587) q[2];
sx q[2];
rz(-3.0465872) q[2];
rz(0.57560086) q[3];
sx q[3];
rz(-2.3592981) q[3];
sx q[3];
rz(1.4466205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59548241) q[0];
sx q[0];
rz(-0.70850104) q[0];
sx q[0];
rz(0.27221671) q[0];
rz(-3.0369924) q[1];
sx q[1];
rz(-2.101254) q[1];
sx q[1];
rz(-1.4629755) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.187785) q[0];
sx q[0];
rz(-1.6208795) q[0];
sx q[0];
rz(-0.84651504) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5065423) q[2];
sx q[2];
rz(-2.3550866) q[2];
sx q[2];
rz(2.3960339) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.21468563) q[1];
sx q[1];
rz(-0.72539293) q[1];
sx q[1];
rz(-2.7100485) q[1];
rz(-pi) q[2];
rz(2.9862254) q[3];
sx q[3];
rz(-0.37305635) q[3];
sx q[3];
rz(-3.1257926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0766803) q[2];
sx q[2];
rz(-1.400759) q[2];
sx q[2];
rz(2.7317375) q[2];
rz(-0.05154933) q[3];
sx q[3];
rz(-0.99286538) q[3];
sx q[3];
rz(0.73369098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9728397) q[0];
sx q[0];
rz(-0.77744716) q[0];
sx q[0];
rz(-0.69946104) q[0];
rz(0.93060023) q[1];
sx q[1];
rz(-1.365463) q[1];
sx q[1];
rz(2.0420989) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4874026) q[0];
sx q[0];
rz(-0.64399566) q[0];
sx q[0];
rz(-1.0985159) q[0];
rz(-2.900057) q[2];
sx q[2];
rz(-2.1061828) q[2];
sx q[2];
rz(0.33171847) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.857497) q[1];
sx q[1];
rz(-1.2624143) q[1];
sx q[1];
rz(-2.8055138) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52249281) q[3];
sx q[3];
rz(-2.429649) q[3];
sx q[3];
rz(-2.5639736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3804271) q[2];
sx q[2];
rz(-2.0695504) q[2];
sx q[2];
rz(2.2011444) q[2];
rz(2.579651) q[3];
sx q[3];
rz(-2.1886531) q[3];
sx q[3];
rz(-0.57653069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024427323) q[0];
sx q[0];
rz(-3.0553525) q[0];
sx q[0];
rz(1.0166136) q[0];
rz(0.92574614) q[1];
sx q[1];
rz(-0.75030202) q[1];
sx q[1];
rz(-3.064916) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0177855) q[0];
sx q[0];
rz(-1.9636256) q[0];
sx q[0];
rz(2.9621302) q[0];
x q[1];
rz(-0.87575298) q[2];
sx q[2];
rz(-0.31538439) q[2];
sx q[2];
rz(-1.8502667) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8342585) q[1];
sx q[1];
rz(-1.3594207) q[1];
sx q[1];
rz(-2.1366875) q[1];
rz(1.881095) q[3];
sx q[3];
rz(-2.3633133) q[3];
sx q[3];
rz(-2.8097514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.76806796) q[2];
sx q[2];
rz(-2.3536451) q[2];
sx q[2];
rz(-0.12316556) q[2];
rz(-2.4672616) q[3];
sx q[3];
rz(-0.47334039) q[3];
sx q[3];
rz(-2.3136852) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3097565) q[0];
sx q[0];
rz(-0.18853822) q[0];
sx q[0];
rz(-2.6605666) q[0];
rz(-2.9006529) q[1];
sx q[1];
rz(-0.53673458) q[1];
sx q[1];
rz(3.0066838) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1447159) q[0];
sx q[0];
rz(-2.2149388) q[0];
sx q[0];
rz(0.5818602) q[0];
x q[1];
rz(2.606484) q[2];
sx q[2];
rz(-1.6596748) q[2];
sx q[2];
rz(-0.8366226) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72249216) q[1];
sx q[1];
rz(-2.3416391) q[1];
sx q[1];
rz(-1.1833619) q[1];
rz(-2.6916299) q[3];
sx q[3];
rz(-1.5057766) q[3];
sx q[3];
rz(-1.7700333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2065108) q[2];
sx q[2];
rz(-0.93588459) q[2];
sx q[2];
rz(-1.3235462) q[2];
rz(0.27213085) q[3];
sx q[3];
rz(-2.2836253) q[3];
sx q[3];
rz(2.9959397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018933522) q[0];
sx q[0];
rz(-0.7974565) q[0];
sx q[0];
rz(-0.65852273) q[0];
rz(-1.5497442) q[1];
sx q[1];
rz(-1.0127944) q[1];
sx q[1];
rz(-0.50643593) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50473266) q[0];
sx q[0];
rz(-1.2018645) q[0];
sx q[0];
rz(1.641666) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1108702) q[2];
sx q[2];
rz(-1.0892158) q[2];
sx q[2];
rz(0.55847634) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5665019) q[1];
sx q[1];
rz(-3.021214) q[1];
sx q[1];
rz(-0.68415227) q[1];
x q[2];
rz(-3.0783404) q[3];
sx q[3];
rz(-2.5613385) q[3];
sx q[3];
rz(1.3903416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8062313) q[2];
sx q[2];
rz(-0.74593097) q[2];
sx q[2];
rz(-1.8719505) q[2];
rz(1.7757724) q[3];
sx q[3];
rz(-0.013805496) q[3];
sx q[3];
rz(2.5891916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7954471) q[0];
sx q[0];
rz(-0.2427635) q[0];
sx q[0];
rz(-0.41918293) q[0];
rz(-0.090713352) q[1];
sx q[1];
rz(-2.2234629) q[1];
sx q[1];
rz(-2.1144287) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61495435) q[0];
sx q[0];
rz(-1.1246343) q[0];
sx q[0];
rz(-0.63132186) q[0];
rz(-0.33454169) q[2];
sx q[2];
rz(-3.1174264) q[2];
sx q[2];
rz(-2.6466359) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.60909437) q[1];
sx q[1];
rz(-0.23819345) q[1];
sx q[1];
rz(1.1042117) q[1];
x q[2];
rz(0.88524359) q[3];
sx q[3];
rz(-1.5285057) q[3];
sx q[3];
rz(-1.2219747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.98296982) q[2];
sx q[2];
rz(-2.136844) q[2];
sx q[2];
rz(2.1515004) q[2];
rz(0.31791911) q[3];
sx q[3];
rz(-0.81219321) q[3];
sx q[3];
rz(-3.1032491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7541499) q[0];
sx q[0];
rz(-2.4328572) q[0];
sx q[0];
rz(2.5842174) q[0];
rz(-0.36224657) q[1];
sx q[1];
rz(-1.7049512) q[1];
sx q[1];
rz(2.974496) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8734735) q[0];
sx q[0];
rz(-2.7303006) q[0];
sx q[0];
rz(2.4231572) q[0];
x q[1];
rz(1.3873151) q[2];
sx q[2];
rz(-2.263732) q[2];
sx q[2];
rz(-0.82833457) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5846338) q[1];
sx q[1];
rz(-1.2642197) q[1];
sx q[1];
rz(0.68309288) q[1];
rz(-pi) q[2];
rz(0.07463712) q[3];
sx q[3];
rz(-0.88866808) q[3];
sx q[3];
rz(2.9934237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.96357137) q[2];
sx q[2];
rz(-0.19266291) q[2];
sx q[2];
rz(-2.7131405) q[2];
rz(2.4109449) q[3];
sx q[3];
rz(-2.0615536) q[3];
sx q[3];
rz(2.54125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927602) q[0];
sx q[0];
rz(-1.4880143) q[0];
sx q[0];
rz(-2.0220508) q[0];
rz(-0.66850942) q[1];
sx q[1];
rz(-2.5709277) q[1];
sx q[1];
rz(-0.25984919) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90454532) q[0];
sx q[0];
rz(-1.0759085) q[0];
sx q[0];
rz(-0.2965692) q[0];
rz(0.38358263) q[2];
sx q[2];
rz(-2.0077939) q[2];
sx q[2];
rz(1.1326157) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.47626859) q[1];
sx q[1];
rz(-0.55108738) q[1];
sx q[1];
rz(0.70161201) q[1];
x q[2];
rz(-1.7612475) q[3];
sx q[3];
rz(-0.69882876) q[3];
sx q[3];
rz(0.33635283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4708289) q[2];
sx q[2];
rz(-1.908611) q[2];
sx q[2];
rz(1.0864351) q[2];
rz(-0.49232617) q[3];
sx q[3];
rz(-2.2958675) q[3];
sx q[3];
rz(0.72211784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54871854) q[0];
sx q[0];
rz(-1.508779) q[0];
sx q[0];
rz(-1.4887703) q[0];
rz(1.8863574) q[1];
sx q[1];
rz(-1.3175169) q[1];
sx q[1];
rz(-3.0138737) q[1];
rz(0.095070953) q[2];
sx q[2];
rz(-1.9683206) q[2];
sx q[2];
rz(-1.5299601) q[2];
rz(1.6519573) q[3];
sx q[3];
rz(-2.4774144) q[3];
sx q[3];
rz(-1.7437205) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
