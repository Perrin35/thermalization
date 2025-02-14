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
rz(0.59347403) q[1];
sx q[1];
rz(3.7945336) q[1];
sx q[1];
rz(8.3795587) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0526177) q[0];
sx q[0];
rz(-2.0988905) q[0];
sx q[0];
rz(-1.0333956) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.521131) q[2];
sx q[2];
rz(-0.6383183) q[2];
sx q[2];
rz(0.82439232) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7683623) q[1];
sx q[1];
rz(-1.4700031) q[1];
sx q[1];
rz(-1.412315) q[1];
rz(1.0020489) q[3];
sx q[3];
rz(-1.4706597) q[3];
sx q[3];
rz(1.1796234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.97592252) q[2];
sx q[2];
rz(-2.133635) q[2];
sx q[2];
rz(2.9284076) q[2];
rz(-0.91607696) q[3];
sx q[3];
rz(-2.7824184) q[3];
sx q[3];
rz(-0.38755125) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2857472) q[0];
sx q[0];
rz(-2.2540932) q[0];
sx q[0];
rz(-0.47072738) q[0];
rz(-2.0384906) q[1];
sx q[1];
rz(-2.6328937) q[1];
sx q[1];
rz(-0.57977605) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8743917) q[0];
sx q[0];
rz(-1.6306595) q[0];
sx q[0];
rz(0.52459985) q[0];
x q[1];
rz(-2.6227399) q[2];
sx q[2];
rz(-0.60728549) q[2];
sx q[2];
rz(2.3407206) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4976385) q[1];
sx q[1];
rz(-2.4259287) q[1];
sx q[1];
rz(2.1791451) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1768119) q[3];
sx q[3];
rz(-2.7005115) q[3];
sx q[3];
rz(-1.4853322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1595999) q[2];
sx q[2];
rz(-2.2856568) q[2];
sx q[2];
rz(3.0465872) q[2];
rz(0.57560086) q[3];
sx q[3];
rz(-2.3592981) q[3];
sx q[3];
rz(1.4466205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59548241) q[0];
sx q[0];
rz(-2.4330916) q[0];
sx q[0];
rz(-2.8693759) q[0];
rz(3.0369924) q[1];
sx q[1];
rz(-2.101254) q[1];
sx q[1];
rz(1.4629755) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42727106) q[0];
sx q[0];
rz(-0.84762379) q[0];
sx q[0];
rz(-0.066825213) q[0];
x q[1];
rz(3.077329) q[2];
sx q[2];
rz(-0.78635629) q[2];
sx q[2];
rz(-0.83640316) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.453664) q[1];
sx q[1];
rz(-1.2896104) q[1];
sx q[1];
rz(2.4635876) q[1];
rz(-pi) q[2];
rz(2.7726309) q[3];
sx q[3];
rz(-1.6272244) q[3];
sx q[3];
rz(1.4101613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0766803) q[2];
sx q[2];
rz(-1.400759) q[2];
sx q[2];
rz(-2.7317375) q[2];
rz(0.05154933) q[3];
sx q[3];
rz(-0.99286538) q[3];
sx q[3];
rz(-0.73369098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16875295) q[0];
sx q[0];
rz(-2.3641455) q[0];
sx q[0];
rz(2.4421316) q[0];
rz(2.2109924) q[1];
sx q[1];
rz(-1.365463) q[1];
sx q[1];
rz(-2.0420989) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4874026) q[0];
sx q[0];
rz(-2.497597) q[0];
sx q[0];
rz(-1.0985159) q[0];
x q[1];
rz(2.1191988) q[2];
sx q[2];
rz(-1.363593) q[2];
sx q[2];
rz(-1.1140547) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2840957) q[1];
sx q[1];
rz(-1.2624143) q[1];
sx q[1];
rz(-2.8055138) q[1];
x q[2];
rz(2.6190998) q[3];
sx q[3];
rz(-0.71194369) q[3];
sx q[3];
rz(2.5639736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3804271) q[2];
sx q[2];
rz(-1.0720422) q[2];
sx q[2];
rz(2.2011444) q[2];
rz(-2.579651) q[3];
sx q[3];
rz(-0.95293957) q[3];
sx q[3];
rz(-0.57653069) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024427323) q[0];
sx q[0];
rz(-0.086240135) q[0];
sx q[0];
rz(2.1249791) q[0];
rz(2.2158465) q[1];
sx q[1];
rz(-0.75030202) q[1];
sx q[1];
rz(-0.076676682) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12380713) q[0];
sx q[0];
rz(-1.9636256) q[0];
sx q[0];
rz(2.9621302) q[0];
x q[1];
rz(-0.20598866) q[2];
sx q[2];
rz(-1.3302557) q[2];
sx q[2];
rz(0.57130106) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3073342) q[1];
sx q[1];
rz(-1.7821719) q[1];
sx q[1];
rz(-2.1366875) q[1];
rz(-pi) q[2];
rz(-0.29239817) q[3];
sx q[3];
rz(-0.83857036) q[3];
sx q[3];
rz(-0.75493073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3735247) q[2];
sx q[2];
rz(-2.3536451) q[2];
sx q[2];
rz(-3.0184271) q[2];
rz(-0.67433107) q[3];
sx q[3];
rz(-0.47334039) q[3];
sx q[3];
rz(2.3136852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3097565) q[0];
sx q[0];
rz(-2.9530544) q[0];
sx q[0];
rz(0.48102608) q[0];
rz(2.9006529) q[1];
sx q[1];
rz(-0.53673458) q[1];
sx q[1];
rz(0.13490881) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8277934) q[0];
sx q[0];
rz(-0.83923566) q[0];
sx q[0];
rz(-2.2025432) q[0];
rz(-0.1730072) q[2];
sx q[2];
rz(-2.5998625) q[2];
sx q[2];
rz(0.8828021) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4191005) q[1];
sx q[1];
rz(-0.79995352) q[1];
sx q[1];
rz(-1.9582307) q[1];
rz(-pi) q[2];
rz(-1.4986131) q[3];
sx q[3];
rz(-2.0197387) q[3];
sx q[3];
rz(-0.23060966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2065108) q[2];
sx q[2];
rz(-2.2057081) q[2];
sx q[2];
rz(-1.8180465) q[2];
rz(0.27213085) q[3];
sx q[3];
rz(-0.85796732) q[3];
sx q[3];
rz(-2.9959397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.018933522) q[0];
sx q[0];
rz(-0.7974565) q[0];
sx q[0];
rz(0.65852273) q[0];
rz(1.5918484) q[1];
sx q[1];
rz(-2.1287983) q[1];
sx q[1];
rz(0.50643593) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50473266) q[0];
sx q[0];
rz(-1.9397282) q[0];
sx q[0];
rz(-1.641666) q[0];
x q[1];
rz(2.0525706) q[2];
sx q[2];
rz(-1.5980236) q[2];
sx q[2];
rz(-2.1435062) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5665019) q[1];
sx q[1];
rz(-3.021214) q[1];
sx q[1];
rz(2.4574404) q[1];
x q[2];
rz(1.6122088) q[3];
sx q[3];
rz(-0.99185252) q[3];
sx q[3];
rz(-1.675663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8062313) q[2];
sx q[2];
rz(-2.3956617) q[2];
sx q[2];
rz(1.2696421) q[2];
rz(-1.3658203) q[3];
sx q[3];
rz(-0.013805496) q[3];
sx q[3];
rz(2.5891916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7954471) q[0];
sx q[0];
rz(-2.8988291) q[0];
sx q[0];
rz(2.7224097) q[0];
rz(-3.0508793) q[1];
sx q[1];
rz(-0.9181298) q[1];
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
rz(2.5102708) q[0];
rz(-pi) q[1];
rz(2.807051) q[2];
sx q[2];
rz(-3.1174264) q[2];
sx q[2];
rz(0.4949567) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5324983) q[1];
sx q[1];
rz(-2.9033992) q[1];
sx q[1];
rz(2.0373809) q[1];
rz(-0.88524359) q[3];
sx q[3];
rz(-1.613087) q[3];
sx q[3];
rz(-1.2219747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98296982) q[2];
sx q[2];
rz(-2.136844) q[2];
sx q[2];
rz(-2.1515004) q[2];
rz(0.31791911) q[3];
sx q[3];
rz(-0.81219321) q[3];
sx q[3];
rz(-3.1032491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38744277) q[0];
sx q[0];
rz(-2.4328572) q[0];
sx q[0];
rz(-0.55737525) q[0];
rz(0.36224657) q[1];
sx q[1];
rz(-1.7049512) q[1];
sx q[1];
rz(-2.974496) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8734735) q[0];
sx q[0];
rz(-2.7303006) q[0];
sx q[0];
rz(-2.4231572) q[0];
rz(-pi) q[1];
rz(-1.3873151) q[2];
sx q[2];
rz(-0.87786061) q[2];
sx q[2];
rz(-0.82833457) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9145979) q[1];
sx q[1];
rz(-0.92508864) q[1];
sx q[1];
rz(1.1832994) q[1];
x q[2];
rz(-1.4792419) q[3];
sx q[3];
rz(-2.4560438) q[3];
sx q[3];
rz(2.8753672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.96357137) q[2];
sx q[2];
rz(-0.19266291) q[2];
sx q[2];
rz(-0.42845217) q[2];
rz(2.4109449) q[3];
sx q[3];
rz(-1.080039) q[3];
sx q[3];
rz(-2.54125) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6927602) q[0];
sx q[0];
rz(-1.4880143) q[0];
sx q[0];
rz(1.1195419) q[0];
rz(0.66850942) q[1];
sx q[1];
rz(-2.5709277) q[1];
sx q[1];
rz(0.25984919) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90454532) q[0];
sx q[0];
rz(-1.0759085) q[0];
sx q[0];
rz(-2.8450235) q[0];
rz(-pi) q[1];
rz(-0.89533676) q[2];
sx q[2];
rz(-0.5731715) q[2];
sx q[2];
rz(1.8942539) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7185134) q[1];
sx q[1];
rz(-1.2260409) q[1];
sx q[1];
rz(-0.43889795) q[1];
rz(0.15774653) q[3];
sx q[3];
rz(-0.88706568) q[3];
sx q[3];
rz(0.58303787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4708289) q[2];
sx q[2];
rz(-1.908611) q[2];
sx q[2];
rz(-1.0864351) q[2];
rz(2.6492665) q[3];
sx q[3];
rz(-0.84572518) q[3];
sx q[3];
rz(2.4194748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
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
rz(1.9699388) q[2];
sx q[2];
rz(-1.4831586) q[2];
sx q[2];
rz(-3.0638564) q[2];
rz(-1.4896354) q[3];
sx q[3];
rz(-2.4774144) q[3];
sx q[3];
rz(-1.7437205) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
