OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.90091997) q[0];
sx q[0];
rz(-0.14578851) q[0];
sx q[0];
rz(-1.8426275) q[0];
rz(2.9453912) q[1];
sx q[1];
rz(-1.8008404) q[1];
sx q[1];
rz(0.10032108) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8226877) q[0];
sx q[0];
rz(-1.5189369) q[0];
sx q[0];
rz(2.479631) q[0];
rz(-pi) q[1];
rz(2.6703628) q[2];
sx q[2];
rz(-1.2696075) q[2];
sx q[2];
rz(-2.9848841) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8954017) q[1];
sx q[1];
rz(-2.7936087) q[1];
sx q[1];
rz(2.4892508) q[1];
rz(-pi) q[2];
rz(-2.6907608) q[3];
sx q[3];
rz(-0.63967645) q[3];
sx q[3];
rz(1.0420639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4002865) q[2];
sx q[2];
rz(-2.9897959) q[2];
sx q[2];
rz(-2.3347704) q[2];
rz(0.72871366) q[3];
sx q[3];
rz(-2.3829134) q[3];
sx q[3];
rz(-1.9369283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.62419409) q[0];
sx q[0];
rz(-0.26573467) q[0];
sx q[0];
rz(-2.8345795) q[0];
rz(-1.2767876) q[1];
sx q[1];
rz(-2.0025608) q[1];
sx q[1];
rz(-2.0397287) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66070304) q[0];
sx q[0];
rz(-2.1212949) q[0];
sx q[0];
rz(-2.9376229) q[0];
rz(-1.6346143) q[2];
sx q[2];
rz(-1.7498657) q[2];
sx q[2];
rz(-0.049953559) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9872322) q[1];
sx q[1];
rz(-1.7061966) q[1];
sx q[1];
rz(1.3920648) q[1];
rz(-pi) q[2];
rz(-2.3526221) q[3];
sx q[3];
rz(-1.0010825) q[3];
sx q[3];
rz(2.8017442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.028194204) q[2];
sx q[2];
rz(-0.58176175) q[2];
sx q[2];
rz(0.096435189) q[2];
rz(0.15245572) q[3];
sx q[3];
rz(-1.5072482) q[3];
sx q[3];
rz(-2.4436387) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089712791) q[0];
sx q[0];
rz(-2.0413601) q[0];
sx q[0];
rz(-0.40019792) q[0];
rz(-1.5125037) q[1];
sx q[1];
rz(-0.17833231) q[1];
sx q[1];
rz(0.2581183) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011964762) q[0];
sx q[0];
rz(-1.3370378) q[0];
sx q[0];
rz(0.14272228) q[0];
rz(0.67419184) q[2];
sx q[2];
rz(-0.12756824) q[2];
sx q[2];
rz(1.3889165) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5096881) q[1];
sx q[1];
rz(-0.013876112) q[1];
sx q[1];
rz(-0.3664151) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24127023) q[3];
sx q[3];
rz(-1.4033969) q[3];
sx q[3];
rz(-0.77360717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1199946) q[2];
sx q[2];
rz(-1.1979411) q[2];
sx q[2];
rz(3.1094587) q[2];
rz(2.8739127) q[3];
sx q[3];
rz(-1.5175502) q[3];
sx q[3];
rz(-1.5599498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5716008) q[0];
sx q[0];
rz(-1.6331693) q[0];
sx q[0];
rz(0.084240325) q[0];
rz(0.035942297) q[1];
sx q[1];
rz(-3.1075931) q[1];
sx q[1];
rz(2.8004004) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55732176) q[0];
sx q[0];
rz(-0.94573089) q[0];
sx q[0];
rz(1.2465501) q[0];
rz(-pi) q[1];
rz(0.3772328) q[2];
sx q[2];
rz(-0.62059852) q[2];
sx q[2];
rz(-1.1334186) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.90662876) q[1];
sx q[1];
rz(-0.90032691) q[1];
sx q[1];
rz(-0.57256289) q[1];
rz(-0.59935948) q[3];
sx q[3];
rz(-1.9512358) q[3];
sx q[3];
rz(-2.3886556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.42698947) q[2];
sx q[2];
rz(-2.0689071) q[2];
sx q[2];
rz(1.535447) q[2];
rz(2.8793907) q[3];
sx q[3];
rz(-1.5235135) q[3];
sx q[3];
rz(-2.1867627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3839805) q[0];
sx q[0];
rz(-2.7181427) q[0];
sx q[0];
rz(2.0641932) q[0];
rz(0.39240882) q[1];
sx q[1];
rz(-0.078374021) q[1];
sx q[1];
rz(2.1108625) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01973132) q[0];
sx q[0];
rz(-1.9872338) q[0];
sx q[0];
rz(-0.58953489) q[0];
x q[1];
rz(1.2472146) q[2];
sx q[2];
rz(-2.6204118) q[2];
sx q[2];
rz(-0.22400907) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5009119) q[1];
sx q[1];
rz(-1.3488773) q[1];
sx q[1];
rz(0.16900105) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34319539) q[3];
sx q[3];
rz(-2.0569909) q[3];
sx q[3];
rz(-2.9064158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82421676) q[2];
sx q[2];
rz(-0.65418303) q[2];
sx q[2];
rz(-2.3023494) q[2];
rz(-0.9134891) q[3];
sx q[3];
rz(-1.8279816) q[3];
sx q[3];
rz(2.998108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4588673) q[0];
sx q[0];
rz(-0.23696466) q[0];
sx q[0];
rz(1.6656026) q[0];
rz(2.7492375) q[1];
sx q[1];
rz(-1.0959492) q[1];
sx q[1];
rz(2.5501693) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2570193) q[0];
sx q[0];
rz(-1.9396175) q[0];
sx q[0];
rz(0.9414282) q[0];
rz(-2.8663581) q[2];
sx q[2];
rz(-1.7448336) q[2];
sx q[2];
rz(1.7013719) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.082669584) q[1];
sx q[1];
rz(-1.5048522) q[1];
sx q[1];
rz(2.0772499) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7631432) q[3];
sx q[3];
rz(-0.90988084) q[3];
sx q[3];
rz(1.701041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8412987) q[2];
sx q[2];
rz(-0.66606194) q[2];
sx q[2];
rz(0.57233125) q[2];
rz(0.22219292) q[3];
sx q[3];
rz(-2.7100345) q[3];
sx q[3];
rz(0.64479327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38729024) q[0];
sx q[0];
rz(-0.13986762) q[0];
sx q[0];
rz(0.40827665) q[0];
rz(-0.73221842) q[1];
sx q[1];
rz(-0.1258985) q[1];
sx q[1];
rz(-2.8439723) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10935303) q[0];
sx q[0];
rz(-1.9178357) q[0];
sx q[0];
rz(-0.45022398) q[0];
rz(-3.1213644) q[2];
sx q[2];
rz(-2.0654404) q[2];
sx q[2];
rz(1.1927562) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0875594) q[1];
sx q[1];
rz(-0.39759025) q[1];
sx q[1];
rz(-1.041574) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70370378) q[3];
sx q[3];
rz(-1.6856442) q[3];
sx q[3];
rz(2.5435257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.97062651) q[2];
sx q[2];
rz(-1.2622702) q[2];
sx q[2];
rz(-2.4453898) q[2];
rz(-0.89037406) q[3];
sx q[3];
rz(-1.1768769) q[3];
sx q[3];
rz(1.4674998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9625229) q[0];
sx q[0];
rz(-3.1130377) q[0];
sx q[0];
rz(0.21275511) q[0];
rz(-2.6720324) q[1];
sx q[1];
rz(-2.1766365) q[1];
sx q[1];
rz(-0.75417095) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.629166) q[0];
sx q[0];
rz(-2.7943724) q[0];
sx q[0];
rz(0.99033333) q[0];
x q[1];
rz(-2.8163359) q[2];
sx q[2];
rz(-1.7647247) q[2];
sx q[2];
rz(3.0162899) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.27867815) q[1];
sx q[1];
rz(-2.2154741) q[1];
sx q[1];
rz(1.7134604) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2175104) q[3];
sx q[3];
rz(-1.0093186) q[3];
sx q[3];
rz(-1.9677066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0065877) q[2];
sx q[2];
rz(-2.2392515) q[2];
sx q[2];
rz(2.3593486) q[2];
rz(-1.4462645) q[3];
sx q[3];
rz(-2.5952314) q[3];
sx q[3];
rz(-2.7915891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0924454) q[0];
sx q[0];
rz(-2.6706084) q[0];
sx q[0];
rz(2.1771722) q[0];
rz(1.8655221) q[1];
sx q[1];
rz(-1.7274906) q[1];
sx q[1];
rz(-1.6395578) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1196331) q[0];
sx q[0];
rz(-0.30888452) q[0];
sx q[0];
rz(0.85794411) q[0];
x q[1];
rz(2.4110498) q[2];
sx q[2];
rz(-1.2337832) q[2];
sx q[2];
rz(2.4180129) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6148541) q[1];
sx q[1];
rz(-2.6236218) q[1];
sx q[1];
rz(1.5244085) q[1];
rz(1.9988868) q[3];
sx q[3];
rz(-2.9799298) q[3];
sx q[3];
rz(2.3695994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8883349) q[2];
sx q[2];
rz(-1.2829245) q[2];
sx q[2];
rz(-2.1962732) q[2];
rz(-2.3156598) q[3];
sx q[3];
rz(-1.632894) q[3];
sx q[3];
rz(-0.53502423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.076040529) q[0];
sx q[0];
rz(-1.9726418) q[0];
sx q[0];
rz(0.8031351) q[0];
rz(1.5777292) q[1];
sx q[1];
rz(-1.6604661) q[1];
sx q[1];
rz(0.28958431) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6824376) q[0];
sx q[0];
rz(-1.4652325) q[0];
sx q[0];
rz(-3.1292874) q[0];
x q[1];
rz(2.824823) q[2];
sx q[2];
rz(-1.8837591) q[2];
sx q[2];
rz(-0.53887689) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3037422) q[1];
sx q[1];
rz(-1.9281862) q[1];
sx q[1];
rz(-0.2405432) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8554162) q[3];
sx q[3];
rz(-1.2499785) q[3];
sx q[3];
rz(2.8367354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3760066) q[2];
sx q[2];
rz(-0.12066081) q[2];
sx q[2];
rz(2.181459) q[2];
rz(0.58297408) q[3];
sx q[3];
rz(-0.65810242) q[3];
sx q[3];
rz(0.8647024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.690602) q[0];
sx q[0];
rz(-1.4324181) q[0];
sx q[0];
rz(1.6313534) q[0];
rz(0.040738978) q[1];
sx q[1];
rz(-0.67650411) q[1];
sx q[1];
rz(0.13112851) q[1];
rz(-1.1004943) q[2];
sx q[2];
rz(-1.8467156) q[2];
sx q[2];
rz(2.4450532) q[2];
rz(-1.4996281) q[3];
sx q[3];
rz(-1.5200079) q[3];
sx q[3];
rz(-0.11848371) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
