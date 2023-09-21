OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(-2.5928901) q[0];
sx q[0];
rz(0.8843511) q[0];
rz(-1.7110775) q[1];
sx q[1];
rz(-0.95354748) q[1];
sx q[1];
rz(-1.5024827) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2384773) q[0];
sx q[0];
rz(-1.3912541) q[0];
sx q[0];
rz(1.205501) q[0];
rz(1.5760954) q[2];
sx q[2];
rz(-1.0056579) q[2];
sx q[2];
rz(-1.1331171) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2734387) q[1];
sx q[1];
rz(-1.8167129) q[1];
sx q[1];
rz(1.6383365) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1611657) q[3];
sx q[3];
rz(-0.75244609) q[3];
sx q[3];
rz(-2.9890271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.78645906) q[2];
sx q[2];
rz(-0.81374514) q[2];
sx q[2];
rz(-0.65594977) q[2];
rz(-1.2077228) q[3];
sx q[3];
rz(-1.9658807) q[3];
sx q[3];
rz(0.99457994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8939963) q[0];
sx q[0];
rz(-2.7212454) q[0];
sx q[0];
rz(0.43352747) q[0];
rz(-0.22878376) q[1];
sx q[1];
rz(-2.7119535) q[1];
sx q[1];
rz(-3.1343592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9418966) q[0];
sx q[0];
rz(-1.6640088) q[0];
sx q[0];
rz(-1.9641563) q[0];
rz(-pi) q[1];
rz(2.8480808) q[2];
sx q[2];
rz(-1.8732757) q[2];
sx q[2];
rz(0.40700618) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0821973) q[1];
sx q[1];
rz(-0.42947436) q[1];
sx q[1];
rz(-0.55179623) q[1];
x q[2];
rz(-0.90356566) q[3];
sx q[3];
rz(-1.9575319) q[3];
sx q[3];
rz(-0.29153338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5269512) q[2];
sx q[2];
rz(-0.80792892) q[2];
sx q[2];
rz(-2.4439404) q[2];
rz(-3.0200322) q[3];
sx q[3];
rz(-1.2391042) q[3];
sx q[3];
rz(0.30383032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6737297) q[0];
sx q[0];
rz(-2.2255852) q[0];
sx q[0];
rz(1.3695705) q[0];
rz(-1.9000152) q[1];
sx q[1];
rz(-1.4135655) q[1];
sx q[1];
rz(2.8799768) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81944377) q[0];
sx q[0];
rz(-1.561164) q[0];
sx q[0];
rz(-1.5887512) q[0];
x q[1];
rz(-3.0171266) q[2];
sx q[2];
rz(-2.3111768) q[2];
sx q[2];
rz(-2.894573) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6690327) q[1];
sx q[1];
rz(-1.6918039) q[1];
sx q[1];
rz(2.3585412) q[1];
rz(1.5393799) q[3];
sx q[3];
rz(-1.0580225) q[3];
sx q[3];
rz(-2.153742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4613142) q[2];
sx q[2];
rz(-1.6363982) q[2];
sx q[2];
rz(0.20351163) q[2];
rz(-2.2198548) q[3];
sx q[3];
rz(-1.2676055) q[3];
sx q[3];
rz(-2.8620499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97776425) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(-1.5699566) q[0];
rz(-2.1381901) q[1];
sx q[1];
rz(-1.3137716) q[1];
sx q[1];
rz(1.8932231) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0715863) q[0];
sx q[0];
rz(-1.4741815) q[0];
sx q[0];
rz(-3.0759401) q[0];
x q[1];
rz(-1.3893045) q[2];
sx q[2];
rz(-1.776812) q[2];
sx q[2];
rz(0.69603053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.215938) q[1];
sx q[1];
rz(-2.2207894) q[1];
sx q[1];
rz(1.710612) q[1];
rz(-0.29420935) q[3];
sx q[3];
rz(-0.74511408) q[3];
sx q[3];
rz(1.7780768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1127597) q[2];
sx q[2];
rz(-1.3112105) q[2];
sx q[2];
rz(-0.83703414) q[2];
rz(-1.2083496) q[3];
sx q[3];
rz(-1.2669867) q[3];
sx q[3];
rz(-0.78554955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0060624881) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(2.3663882) q[0];
rz(-2.7397621) q[1];
sx q[1];
rz(-2.1907175) q[1];
sx q[1];
rz(2.2391589) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1989312) q[0];
sx q[0];
rz(-2.5006602) q[0];
sx q[0];
rz(3.065227) q[0];
rz(-pi) q[1];
rz(2.6685153) q[2];
sx q[2];
rz(-0.50191754) q[2];
sx q[2];
rz(0.15854533) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.202995) q[1];
sx q[1];
rz(-2.1122167) q[1];
sx q[1];
rz(-2.2375537) q[1];
rz(1.8875214) q[3];
sx q[3];
rz(-1.5681019) q[3];
sx q[3];
rz(-2.5731784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.19501413) q[2];
sx q[2];
rz(-2.6533551) q[2];
sx q[2];
rz(1.1966594) q[2];
rz(-1.6992016) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(-1.6825914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.8114132) q[0];
sx q[0];
rz(-0.17689642) q[0];
sx q[0];
rz(0.50317558) q[0];
rz(1.4563837) q[1];
sx q[1];
rz(-2.0676985) q[1];
sx q[1];
rz(-2.9398289) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.266298) q[0];
sx q[0];
rz(-1.504997) q[0];
sx q[0];
rz(-1.468303) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6642338) q[2];
sx q[2];
rz(-1.1672033) q[2];
sx q[2];
rz(1.8237643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4060496) q[1];
sx q[1];
rz(-0.65945259) q[1];
sx q[1];
rz(2.4156648) q[1];
rz(2.6334409) q[3];
sx q[3];
rz(-2.3414632) q[3];
sx q[3];
rz(1.2451764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21489828) q[2];
sx q[2];
rz(-1.639067) q[2];
sx q[2];
rz(-0.26724896) q[2];
rz(-0.8231419) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(-0.87583035) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8478407) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(0.057549495) q[0];
rz(1.6607704) q[1];
sx q[1];
rz(-1.7819504) q[1];
sx q[1];
rz(-0.94271359) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1377624) q[0];
sx q[0];
rz(-2.1592405) q[0];
sx q[0];
rz(-2.0501191) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96037453) q[2];
sx q[2];
rz(-0.27563169) q[2];
sx q[2];
rz(-0.20197091) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.47989935) q[1];
sx q[1];
rz(-2.6091895) q[1];
sx q[1];
rz(2.2366621) q[1];
rz(-0.25992486) q[3];
sx q[3];
rz(-2.0689031) q[3];
sx q[3];
rz(0.92344027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0223579) q[2];
sx q[2];
rz(-1.0417577) q[2];
sx q[2];
rz(2.7589202) q[2];
rz(1.0391957) q[3];
sx q[3];
rz(-1.8363876) q[3];
sx q[3];
rz(-2.1634845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4246178) q[0];
sx q[0];
rz(-0.096352339) q[0];
sx q[0];
rz(2.8714645) q[0];
rz(-0.62942901) q[1];
sx q[1];
rz(-2.4286178) q[1];
sx q[1];
rz(0.28392917) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2487508) q[0];
sx q[0];
rz(-0.98309702) q[0];
sx q[0];
rz(-2.6177004) q[0];
rz(-pi) q[1];
rz(2.9291199) q[2];
sx q[2];
rz(-2.5310235) q[2];
sx q[2];
rz(3.0688822) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.24208454) q[1];
sx q[1];
rz(-0.71166066) q[1];
sx q[1];
rz(2.3942024) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1030032) q[3];
sx q[3];
rz(-2.4814197) q[3];
sx q[3];
rz(-0.71036464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5876864) q[2];
sx q[2];
rz(-0.17051414) q[2];
sx q[2];
rz(1.930687) q[2];
rz(2.8816913) q[3];
sx q[3];
rz(-0.62429684) q[3];
sx q[3];
rz(-2.5134145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.07638409) q[0];
sx q[0];
rz(-0.56088352) q[0];
sx q[0];
rz(2.912345) q[0];
rz(0.30300888) q[1];
sx q[1];
rz(-1.7508933) q[1];
sx q[1];
rz(1.680826) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47855908) q[0];
sx q[0];
rz(-1.8543108) q[0];
sx q[0];
rz(1.5525596) q[0];
x q[1];
rz(0.37947189) q[2];
sx q[2];
rz(-1.737829) q[2];
sx q[2];
rz(-0.52969474) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8252392) q[1];
sx q[1];
rz(-1.5931399) q[1];
sx q[1];
rz(-0.52492001) q[1];
x q[2];
rz(1.1838412) q[3];
sx q[3];
rz(-0.82434067) q[3];
sx q[3];
rz(1.6798306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.71172697) q[2];
sx q[2];
rz(-1.2167598) q[2];
sx q[2];
rz(-0.6955859) q[2];
rz(0.43186489) q[3];
sx q[3];
rz(-2.6769107) q[3];
sx q[3];
rz(-0.7152043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.5678976) q[0];
sx q[0];
rz(-1.9298113) q[0];
sx q[0];
rz(0.33690548) q[0];
rz(-2.9341872) q[1];
sx q[1];
rz(-1.0284871) q[1];
sx q[1];
rz(0.38063231) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56363737) q[0];
sx q[0];
rz(-1.6407688) q[0];
sx q[0];
rz(-1.3689343) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62217525) q[2];
sx q[2];
rz(-2.2969349) q[2];
sx q[2];
rz(-3.0058793) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8457348) q[1];
sx q[1];
rz(-1.0746135) q[1];
sx q[1];
rz(2.9198398) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57772824) q[3];
sx q[3];
rz(-1.2168222) q[3];
sx q[3];
rz(-2.3059394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1404861) q[2];
sx q[2];
rz(-1.1662741) q[2];
sx q[2];
rz(2.005119) q[2];
rz(-3.100637) q[3];
sx q[3];
rz(-2.332873) q[3];
sx q[3];
rz(1.827318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3363591) q[0];
sx q[0];
rz(-2.2362066) q[0];
sx q[0];
rz(2.8295828) q[0];
rz(1.0271172) q[1];
sx q[1];
rz(-1.849091) q[1];
sx q[1];
rz(-1.0277933) q[1];
rz(0.74577352) q[2];
sx q[2];
rz(-0.33200982) q[2];
sx q[2];
rz(0.40340323) q[2];
rz(2.0102262) q[3];
sx q[3];
rz(-1.3888748) q[3];
sx q[3];
rz(0.5007762) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
