OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6498123) q[0];
sx q[0];
rz(-0.28591135) q[0];
sx q[0];
rz(0.51529348) q[0];
rz(-1.7973068) q[1];
sx q[1];
rz(-0.15434115) q[1];
sx q[1];
rz(-0.57758346) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5776423) q[0];
sx q[0];
rz(-2.2221774) q[0];
sx q[0];
rz(1.3629859) q[0];
x q[1];
rz(-2.3876786) q[2];
sx q[2];
rz(-2.4587817) q[2];
sx q[2];
rz(-0.66893259) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2204199) q[1];
sx q[1];
rz(-2.6487659) q[1];
sx q[1];
rz(-0.9763896) q[1];
rz(2.4926315) q[3];
sx q[3];
rz(-2.0994086) q[3];
sx q[3];
rz(1.987207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43705964) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(2.4543767) q[2];
rz(-1.0152738) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(3.0190873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9706443) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(-1.2600391) q[0];
rz(-1.0062224) q[1];
sx q[1];
rz(-2.1496014) q[1];
sx q[1];
rz(-2.2959183) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2226919) q[0];
sx q[0];
rz(-2.0881623) q[0];
sx q[0];
rz(-1.067724) q[0];
x q[1];
rz(-0.77180736) q[2];
sx q[2];
rz(-0.1856005) q[2];
sx q[2];
rz(-1.1112569) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0963124) q[1];
sx q[1];
rz(-1.7573866) q[1];
sx q[1];
rz(0.50218302) q[1];
x q[2];
rz(-2.2250697) q[3];
sx q[3];
rz(-0.73361165) q[3];
sx q[3];
rz(-0.77378002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3530897) q[2];
sx q[2];
rz(-0.22558364) q[2];
sx q[2];
rz(2.6611924) q[2];
rz(1.3530312) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(-1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1903494) q[0];
sx q[0];
rz(-0.23878637) q[0];
sx q[0];
rz(-2.3685266) q[0];
rz(3.0103325) q[1];
sx q[1];
rz(-1.8570329) q[1];
sx q[1];
rz(-1.0864331) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4341136) q[0];
sx q[0];
rz(-1.1093603) q[0];
sx q[0];
rz(3.1397318) q[0];
rz(-pi) q[1];
rz(-1.7227206) q[2];
sx q[2];
rz(-0.65290367) q[2];
sx q[2];
rz(0.34130794) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4194581) q[1];
sx q[1];
rz(-0.28885435) q[1];
sx q[1];
rz(0.018317776) q[1];
rz(-pi) q[2];
rz(-1.1826035) q[3];
sx q[3];
rz(-0.43366323) q[3];
sx q[3];
rz(-0.38503034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0125668) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(1.770299) q[2];
rz(0.38315547) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(-2.3390521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6894158) q[0];
sx q[0];
rz(-1.8912264) q[0];
sx q[0];
rz(3.0932328) q[0];
rz(-2.9776749) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(-1.4455459) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7955129) q[0];
sx q[0];
rz(-0.68329408) q[0];
sx q[0];
rz(0.28738316) q[0];
rz(-pi) q[1];
rz(-1.4807329) q[2];
sx q[2];
rz(-1.3552595) q[2];
sx q[2];
rz(-0.71722523) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1452892) q[1];
sx q[1];
rz(-0.22876303) q[1];
sx q[1];
rz(0.28782515) q[1];
rz(-1.0501782) q[3];
sx q[3];
rz(-1.8225267) q[3];
sx q[3];
rz(0.43581918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0754898) q[2];
sx q[2];
rz(-1.3572071) q[2];
sx q[2];
rz(-1.0243105) q[2];
rz(-1.5284437) q[3];
sx q[3];
rz(-1.5214835) q[3];
sx q[3];
rz(-2.9083692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7016474) q[0];
sx q[0];
rz(-2.3174536) q[0];
sx q[0];
rz(1.2874999) q[0];
rz(-0.31907407) q[1];
sx q[1];
rz(-1.5998452) q[1];
sx q[1];
rz(-2.2873926) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4860977) q[0];
sx q[0];
rz(-2.4276519) q[0];
sx q[0];
rz(-1.5582725) q[0];
x q[1];
rz(-0.87128432) q[2];
sx q[2];
rz(-1.185002) q[2];
sx q[2];
rz(-2.3988349) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2367868) q[1];
sx q[1];
rz(-0.85610897) q[1];
sx q[1];
rz(-2.3805815) q[1];
x q[2];
rz(-0.55027996) q[3];
sx q[3];
rz(-2.0645421) q[3];
sx q[3];
rz(2.9328336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1896818) q[2];
sx q[2];
rz(-0.59331912) q[2];
sx q[2];
rz(0.577315) q[2];
rz(-2.632085) q[3];
sx q[3];
rz(-2.7039492) q[3];
sx q[3];
rz(2.234941) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0034870738) q[0];
sx q[0];
rz(-2.0697937) q[0];
sx q[0];
rz(0.072120897) q[0];
rz(-2.0347319) q[1];
sx q[1];
rz(-2.6289584) q[1];
sx q[1];
rz(0.12621005) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70996767) q[0];
sx q[0];
rz(-2.9266848) q[0];
sx q[0];
rz(-0.14847319) q[0];
rz(0.77504471) q[2];
sx q[2];
rz(-1.2942874) q[2];
sx q[2];
rz(-1.0724049) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5050161) q[1];
sx q[1];
rz(-2.2762183) q[1];
sx q[1];
rz(-0.24800639) q[1];
rz(-1.0522271) q[3];
sx q[3];
rz(-2.7675981) q[3];
sx q[3];
rz(-0.75043375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.90298992) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(-0.77511707) q[2];
rz(0.827968) q[3];
sx q[3];
rz(-0.29100806) q[3];
sx q[3];
rz(-2.0194139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5489952) q[0];
sx q[0];
rz(-2.7153375) q[0];
sx q[0];
rz(-3.0431842) q[0];
rz(-1.1920284) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(-0.55955204) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3768809) q[0];
sx q[0];
rz(-1.6983713) q[0];
sx q[0];
rz(2.5401831) q[0];
rz(-pi) q[1];
rz(-2.9497629) q[2];
sx q[2];
rz(-0.26608135) q[2];
sx q[2];
rz(1.9907469) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2110062) q[1];
sx q[1];
rz(-1.483327) q[1];
sx q[1];
rz(-1.4409815) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3748752) q[3];
sx q[3];
rz(-0.94983263) q[3];
sx q[3];
rz(1.2697112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.45903912) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(-0.34379488) q[2];
rz(2.5750459) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(-2.6678273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.375181) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(-0.73076105) q[0];
rz(2.9991951) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(-2.2699845) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14411892) q[0];
sx q[0];
rz(-1.561957) q[0];
sx q[0];
rz(-3.0661422) q[0];
x q[1];
rz(2.7360104) q[2];
sx q[2];
rz(-1.7103346) q[2];
sx q[2];
rz(1.9987193) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1069113) q[1];
sx q[1];
rz(-2.2195663) q[1];
sx q[1];
rz(-1.1120863) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3685162) q[3];
sx q[3];
rz(-0.61198046) q[3];
sx q[3];
rz(-0.66093854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.72620755) q[2];
sx q[2];
rz(-2.0724847) q[2];
sx q[2];
rz(-2.0020206) q[2];
rz(-1.4987882) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(2.22877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9386439) q[0];
sx q[0];
rz(-1.6904172) q[0];
sx q[0];
rz(-1.2217481) q[0];
rz(2.9755759) q[1];
sx q[1];
rz(-1.8211726) q[1];
sx q[1];
rz(-1.6171914) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4964543) q[0];
sx q[0];
rz(-1.5418515) q[0];
sx q[0];
rz(-1.4883947) q[0];
rz(1.379307) q[2];
sx q[2];
rz(-2.0404437) q[2];
sx q[2];
rz(2.7188403) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.081454885) q[1];
sx q[1];
rz(-2.7739035) q[1];
sx q[1];
rz(1.1170438) q[1];
x q[2];
rz(0.41096656) q[3];
sx q[3];
rz(-2.4604359) q[3];
sx q[3];
rz(-0.70511234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1200072) q[2];
sx q[2];
rz(-1.465613) q[2];
sx q[2];
rz(0.35153708) q[2];
rz(1.0567788) q[3];
sx q[3];
rz(-0.52962279) q[3];
sx q[3];
rz(-0.74469152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4979424) q[0];
sx q[0];
rz(-0.90181667) q[0];
sx q[0];
rz(1.3051916) q[0];
rz(0.38048831) q[1];
sx q[1];
rz(-2.0996129) q[1];
sx q[1];
rz(0.25340733) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6127374) q[0];
sx q[0];
rz(-2.7326072) q[0];
sx q[0];
rz(-1.3091875) q[0];
rz(-1.1240187) q[2];
sx q[2];
rz(-2.7687216) q[2];
sx q[2];
rz(-2.1602221) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0512143) q[1];
sx q[1];
rz(-1.2284632) q[1];
sx q[1];
rz(-1.315209) q[1];
rz(1.9990342) q[3];
sx q[3];
rz(-2.1381452) q[3];
sx q[3];
rz(0.67728562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5499251) q[2];
sx q[2];
rz(-0.88576907) q[2];
sx q[2];
rz(0.5029451) q[2];
rz(-0.89899603) q[3];
sx q[3];
rz(-1.2939724) q[3];
sx q[3];
rz(1.9780654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9713365) q[0];
sx q[0];
rz(-1.6032871) q[0];
sx q[0];
rz(0.26300318) q[0];
rz(2.4304216) q[1];
sx q[1];
rz(-1.0881337) q[1];
sx q[1];
rz(1.7137391) q[1];
rz(-0.74264991) q[2];
sx q[2];
rz(-0.37336083) q[2];
sx q[2];
rz(0.20867418) q[2];
rz(-0.67646277) q[3];
sx q[3];
rz(-2.2623895) q[3];
sx q[3];
rz(-1.1415979) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
