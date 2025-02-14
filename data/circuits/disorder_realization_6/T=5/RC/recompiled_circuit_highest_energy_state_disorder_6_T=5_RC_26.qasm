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
rz(-2.1036086) q[0];
sx q[0];
rz(-1.684364) q[0];
sx q[0];
rz(-1.7888223) q[0];
rz(1.1846722) q[1];
sx q[1];
rz(-0.67263043) q[1];
sx q[1];
rz(-1.6258568) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.420879) q[0];
sx q[0];
rz(-1.7569245) q[0];
sx q[0];
rz(1.3693007) q[0];
x q[1];
rz(2.6625161) q[2];
sx q[2];
rz(-2.276439) q[2];
sx q[2];
rz(0.75032633) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7616854) q[1];
sx q[1];
rz(-2.5406079) q[1];
sx q[1];
rz(0.65324776) q[1];
x q[2];
rz(1.2836841) q[3];
sx q[3];
rz(-1.1653882) q[3];
sx q[3];
rz(1.4472345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1379913) q[2];
sx q[2];
rz(-3.0333952) q[2];
sx q[2];
rz(-0.18599621) q[2];
rz(3.1229535) q[3];
sx q[3];
rz(-1.2942856) q[3];
sx q[3];
rz(-1.0703465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5828534) q[0];
sx q[0];
rz(-2.5753729) q[0];
sx q[0];
rz(-2.8232316) q[0];
rz(-0.63703498) q[1];
sx q[1];
rz(-0.77992264) q[1];
sx q[1];
rz(-0.46924082) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5770188) q[0];
sx q[0];
rz(-1.3718318) q[0];
sx q[0];
rz(1.3964723) q[0];
rz(1.7809243) q[2];
sx q[2];
rz(-0.38345018) q[2];
sx q[2];
rz(-2.0166778) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7982895) q[1];
sx q[1];
rz(-1.8899916) q[1];
sx q[1];
rz(1.809824) q[1];
rz(-2.2607311) q[3];
sx q[3];
rz(-1.5663877) q[3];
sx q[3];
rz(-0.064383163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8239173) q[2];
sx q[2];
rz(-2.9313512) q[2];
sx q[2];
rz(-0.69327411) q[2];
rz(-1.8215826) q[3];
sx q[3];
rz(-1.3258679) q[3];
sx q[3];
rz(-2.0461953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0005118) q[0];
sx q[0];
rz(-0.63758272) q[0];
sx q[0];
rz(-2.6620423) q[0];
rz(-2.6374822) q[1];
sx q[1];
rz(-1.9934374) q[1];
sx q[1];
rz(-3.0832916) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8970015) q[0];
sx q[0];
rz(-2.069733) q[0];
sx q[0];
rz(1.2759491) q[0];
rz(-2.0928755) q[2];
sx q[2];
rz(-2.1914542) q[2];
sx q[2];
rz(-1.716937) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2268277) q[1];
sx q[1];
rz(-1.0221507) q[1];
sx q[1];
rz(0.1424205) q[1];
rz(1.8350868) q[3];
sx q[3];
rz(-1.0796781) q[3];
sx q[3];
rz(2.6437505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38640675) q[2];
sx q[2];
rz(-1.9599954) q[2];
sx q[2];
rz(3.1043261) q[2];
rz(-1.6218328) q[3];
sx q[3];
rz(-0.45791364) q[3];
sx q[3];
rz(2.7601385) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94654361) q[0];
sx q[0];
rz(-2.6730972) q[0];
sx q[0];
rz(-0.066548912) q[0];
rz(-1.5244124) q[1];
sx q[1];
rz(-1.8981372) q[1];
sx q[1];
rz(2.4066511) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83072216) q[0];
sx q[0];
rz(-1.3390216) q[0];
sx q[0];
rz(0.59344296) q[0];
x q[1];
rz(1.8040015) q[2];
sx q[2];
rz(-2.2873023) q[2];
sx q[2];
rz(-0.75762123) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0385602) q[1];
sx q[1];
rz(-1.2913229) q[1];
sx q[1];
rz(0.36906645) q[1];
rz(1.7363343) q[3];
sx q[3];
rz(-2.41417) q[3];
sx q[3];
rz(2.7911428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8526326) q[2];
sx q[2];
rz(-1.5776878) q[2];
sx q[2];
rz(3.0328879) q[2];
rz(-2.2147801) q[3];
sx q[3];
rz(-2.1709397) q[3];
sx q[3];
rz(1.5638117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.0573334) q[0];
sx q[0];
rz(-0.63705343) q[0];
sx q[0];
rz(-2.0507226) q[0];
rz(-0.57251656) q[1];
sx q[1];
rz(-1.9520091) q[1];
sx q[1];
rz(-2.3172839) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.123468) q[0];
sx q[0];
rz(-1.9459007) q[0];
sx q[0];
rz(-2.8897622) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9277497) q[2];
sx q[2];
rz(-2.0111736) q[2];
sx q[2];
rz(1.7945031) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.2356154) q[1];
sx q[1];
rz(-2.4433377) q[1];
sx q[1];
rz(-2.2119616) q[1];
rz(0.82309725) q[3];
sx q[3];
rz(-0.56328008) q[3];
sx q[3];
rz(0.34264229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6606286) q[2];
sx q[2];
rz(-0.97233665) q[2];
sx q[2];
rz(-2.8503897) q[2];
rz(-2.5045942) q[3];
sx q[3];
rz(-1.6773418) q[3];
sx q[3];
rz(1.8891107) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7630735) q[0];
sx q[0];
rz(-2.851649) q[0];
sx q[0];
rz(0.13105233) q[0];
rz(-1.9746926) q[1];
sx q[1];
rz(-0.98375541) q[1];
sx q[1];
rz(0.022620591) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93504923) q[0];
sx q[0];
rz(-0.66842043) q[0];
sx q[0];
rz(-2.5517625) q[0];
rz(2.6220967) q[2];
sx q[2];
rz(-0.24458376) q[2];
sx q[2];
rz(-1.0102538) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.57405942) q[1];
sx q[1];
rz(-1.0242192) q[1];
sx q[1];
rz(-2.778227) q[1];
x q[2];
rz(-1.3778481) q[3];
sx q[3];
rz(-1.7304599) q[3];
sx q[3];
rz(0.99332383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8733946) q[2];
sx q[2];
rz(-2.4384273) q[2];
sx q[2];
rz(-1.084729) q[2];
rz(-1.9966513) q[3];
sx q[3];
rz(-1.8549253) q[3];
sx q[3];
rz(-2.1196608) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5893843) q[0];
sx q[0];
rz(-1.5188058) q[0];
sx q[0];
rz(0.4735221) q[0];
rz(-1.2208968) q[1];
sx q[1];
rz(-0.41243204) q[1];
sx q[1];
rz(-1.268505) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4421688) q[0];
sx q[0];
rz(-1.0528565) q[0];
sx q[0];
rz(0.59364001) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94077295) q[2];
sx q[2];
rz(-2.2029999) q[2];
sx q[2];
rz(-0.8559627) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.699325) q[1];
sx q[1];
rz(-2.0751157) q[1];
sx q[1];
rz(-1.0902576) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0067665) q[3];
sx q[3];
rz(-1.4957168) q[3];
sx q[3];
rz(2.1030542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.97741693) q[2];
sx q[2];
rz(-1.5722534) q[2];
sx q[2];
rz(-0.24641985) q[2];
rz(-0.94270802) q[3];
sx q[3];
rz(-2.0556512) q[3];
sx q[3];
rz(0.76343083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8134269) q[0];
sx q[0];
rz(-1.8712217) q[0];
sx q[0];
rz(3.0810007) q[0];
rz(-1.6391485) q[1];
sx q[1];
rz(-1.363058) q[1];
sx q[1];
rz(-1.5477808) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0133724) q[0];
sx q[0];
rz(-0.0039847535) q[0];
sx q[0];
rz(0.26040034) q[0];
rz(-pi) q[1];
rz(2.5663788) q[2];
sx q[2];
rz(-1.7292495) q[2];
sx q[2];
rz(2.2666933) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.21067737) q[1];
sx q[1];
rz(-2.0077171) q[1];
sx q[1];
rz(2.2411437) q[1];
rz(-pi) q[2];
rz(2.4259858) q[3];
sx q[3];
rz(-1.2168435) q[3];
sx q[3];
rz(-2.7277814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26555201) q[2];
sx q[2];
rz(-2.2206842) q[2];
sx q[2];
rz(-1.6843686) q[2];
rz(-3.0217116) q[3];
sx q[3];
rz(-2.9370152) q[3];
sx q[3];
rz(-2.1129107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2940755) q[0];
sx q[0];
rz(-0.99440044) q[0];
sx q[0];
rz(-2.0062398) q[0];
rz(-3.0382233) q[1];
sx q[1];
rz(-1.7366948) q[1];
sx q[1];
rz(0.9476544) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68610588) q[0];
sx q[0];
rz(-2.4999833) q[0];
sx q[0];
rz(-3.0296586) q[0];
x q[1];
rz(-2.060076) q[2];
sx q[2];
rz(-0.7502509) q[2];
sx q[2];
rz(-2.3231151) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.39571849) q[1];
sx q[1];
rz(-1.304879) q[1];
sx q[1];
rz(-0.87541343) q[1];
rz(-0.89449785) q[3];
sx q[3];
rz(-0.94493659) q[3];
sx q[3];
rz(1.2813527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7171628) q[2];
sx q[2];
rz(-0.62817502) q[2];
sx q[2];
rz(-2.2372645) q[2];
rz(1.2633911) q[3];
sx q[3];
rz(-1.9046013) q[3];
sx q[3];
rz(2.625107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.271027) q[0];
sx q[0];
rz(-2.3834383) q[0];
sx q[0];
rz(0.98536056) q[0];
rz(-0.35161463) q[1];
sx q[1];
rz(-0.96013394) q[1];
sx q[1];
rz(0.34686372) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7642286) q[0];
sx q[0];
rz(-2.1954311) q[0];
sx q[0];
rz(0.2701811) q[0];
rz(-0.11254452) q[2];
sx q[2];
rz(-2.0142609) q[2];
sx q[2];
rz(2.7751768) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.092791768) q[1];
sx q[1];
rz(-1.4818703) q[1];
sx q[1];
rz(-0.6912749) q[1];
rz(-pi) q[2];
rz(-0.44054359) q[3];
sx q[3];
rz(-1.2367354) q[3];
sx q[3];
rz(1.7831217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.88006192) q[2];
sx q[2];
rz(-0.71892771) q[2];
sx q[2];
rz(-0.13492179) q[2];
rz(-3.0123582) q[3];
sx q[3];
rz(-0.40004572) q[3];
sx q[3];
rz(2.1028886) q[3];
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
rz(-pi) q[0];
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
rz(-3.1226817) q[0];
sx q[0];
rz(-1.2553348) q[0];
sx q[0];
rz(-3.0154764) q[0];
rz(-0.46161721) q[1];
sx q[1];
rz(-1.7414265) q[1];
sx q[1];
rz(0.19207676) q[1];
rz(2.6814396) q[2];
sx q[2];
rz(-0.94993409) q[2];
sx q[2];
rz(0.86144907) q[2];
rz(0.57030525) q[3];
sx q[3];
rz(-2.8371596) q[3];
sx q[3];
rz(0.53573487) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
