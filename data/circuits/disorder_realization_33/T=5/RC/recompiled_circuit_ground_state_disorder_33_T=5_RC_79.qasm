OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6991718) q[0];
sx q[0];
rz(-0.81096634) q[0];
sx q[0];
rz(0.45642689) q[0];
rz(2.2189848) q[1];
sx q[1];
rz(2.1906617) q[1];
sx q[1];
rz(9.6074109) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73494324) q[0];
sx q[0];
rz(-1.8514086) q[0];
sx q[0];
rz(0.92241241) q[0];
rz(-2.2721108) q[2];
sx q[2];
rz(-1.1966146) q[2];
sx q[2];
rz(-2.2187198) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6485893) q[1];
sx q[1];
rz(-1.3979288) q[1];
sx q[1];
rz(1.6088617) q[1];
x q[2];
rz(0.62599085) q[3];
sx q[3];
rz(-2.5267793) q[3];
sx q[3];
rz(1.2683319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38723382) q[2];
sx q[2];
rz(-2.0626455) q[2];
sx q[2];
rz(-2.8299502) q[2];
rz(-1.257487) q[3];
sx q[3];
rz(-0.25185549) q[3];
sx q[3];
rz(2.9529412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99609128) q[0];
sx q[0];
rz(-1.6956734) q[0];
sx q[0];
rz(-1.9022994) q[0];
rz(0.39341012) q[1];
sx q[1];
rz(-1.0478123) q[1];
sx q[1];
rz(-2.1107103) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9728972) q[0];
sx q[0];
rz(-2.4542232) q[0];
sx q[0];
rz(2.1363791) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3412881) q[2];
sx q[2];
rz(-1.3131427) q[2];
sx q[2];
rz(2.9837556) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1340902) q[1];
sx q[1];
rz(-0.94408997) q[1];
sx q[1];
rz(-0.43372633) q[1];
rz(-pi) q[2];
rz(0.83194701) q[3];
sx q[3];
rz(-2.0752677) q[3];
sx q[3];
rz(0.17278663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3139412) q[2];
sx q[2];
rz(-0.88259077) q[2];
sx q[2];
rz(-0.35448709) q[2];
rz(0.83141023) q[3];
sx q[3];
rz(-0.76787132) q[3];
sx q[3];
rz(1.6262511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3671234) q[0];
sx q[0];
rz(-2.9614083) q[0];
sx q[0];
rz(-2.6572976) q[0];
rz(-2.2593185) q[1];
sx q[1];
rz(-2.2703998) q[1];
sx q[1];
rz(2.5994515) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4216293) q[0];
sx q[0];
rz(-2.6609592) q[0];
sx q[0];
rz(-0.84893562) q[0];
rz(1.5050921) q[2];
sx q[2];
rz(-1.8915542) q[2];
sx q[2];
rz(1.9707206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.10239359) q[1];
sx q[1];
rz(-0.77966438) q[1];
sx q[1];
rz(1.3179661) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8291953) q[3];
sx q[3];
rz(-0.2451788) q[3];
sx q[3];
rz(2.6441531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.813628) q[2];
sx q[2];
rz(-0.70983228) q[2];
sx q[2];
rz(-1.0343118) q[2];
rz(-1.7799001) q[3];
sx q[3];
rz(-1.1797649) q[3];
sx q[3];
rz(-0.55083197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66083241) q[0];
sx q[0];
rz(-1.7518504) q[0];
sx q[0];
rz(-1.047629) q[0];
rz(-1.818559) q[1];
sx q[1];
rz(-0.58224693) q[1];
sx q[1];
rz(-2.7563162) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40257257) q[0];
sx q[0];
rz(-1.0642576) q[0];
sx q[0];
rz(2.6468524) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0550256) q[2];
sx q[2];
rz(-2.4231152) q[2];
sx q[2];
rz(2.1849146) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.99285728) q[1];
sx q[1];
rz(-2.9210733) q[1];
sx q[1];
rz(-0.33111568) q[1];
rz(-pi) q[2];
rz(1.436211) q[3];
sx q[3];
rz(-2.453605) q[3];
sx q[3];
rz(-2.1768895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72383991) q[2];
sx q[2];
rz(-0.94261348) q[2];
sx q[2];
rz(-2.8670132) q[2];
rz(2.5802021) q[3];
sx q[3];
rz(-2.0761108) q[3];
sx q[3];
rz(1.6339462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0059589) q[0];
sx q[0];
rz(-0.57666403) q[0];
sx q[0];
rz(-2.2744001) q[0];
rz(0.37569702) q[1];
sx q[1];
rz(-2.5521894) q[1];
sx q[1];
rz(-1.2295178) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1389826) q[0];
sx q[0];
rz(-2.0073237) q[0];
sx q[0];
rz(1.3912203) q[0];
x q[1];
rz(-1.0528238) q[2];
sx q[2];
rz(-1.3020294) q[2];
sx q[2];
rz(0.72552272) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2676437) q[1];
sx q[1];
rz(-2.3726844) q[1];
sx q[1];
rz(0.88458367) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2476056) q[3];
sx q[3];
rz(-0.5115307) q[3];
sx q[3];
rz(2.3681896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.025297252) q[2];
sx q[2];
rz(-1.9090434) q[2];
sx q[2];
rz(-3.0787025) q[2];
rz(0.40999117) q[3];
sx q[3];
rz(-2.4153109) q[3];
sx q[3];
rz(-0.15967742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1493688) q[0];
sx q[0];
rz(-0.069644444) q[0];
sx q[0];
rz(-2.3058291) q[0];
rz(-1.958485) q[1];
sx q[1];
rz(-1.2703905) q[1];
sx q[1];
rz(-0.74367181) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0346468) q[0];
sx q[0];
rz(-1.2689044) q[0];
sx q[0];
rz(0.56580122) q[0];
rz(-pi) q[1];
rz(-1.7149107) q[2];
sx q[2];
rz(-1.5169355) q[2];
sx q[2];
rz(0.28648057) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2433155) q[1];
sx q[1];
rz(-0.87496266) q[1];
sx q[1];
rz(-0.23855539) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20009508) q[3];
sx q[3];
rz(-0.26069122) q[3];
sx q[3];
rz(-0.30567769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.30367294) q[2];
sx q[2];
rz(-2.643879) q[2];
sx q[2];
rz(-0.79353235) q[2];
rz(-2.7225336) q[3];
sx q[3];
rz(-1.6520809) q[3];
sx q[3];
rz(0.52217531) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1917052) q[0];
sx q[0];
rz(-0.97923034) q[0];
sx q[0];
rz(1.9784084) q[0];
rz(2.361182) q[1];
sx q[1];
rz(-0.33640948) q[1];
sx q[1];
rz(-1.6132678) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1799729) q[0];
sx q[0];
rz(-0.69937569) q[0];
sx q[0];
rz(1.6716624) q[0];
rz(-0.18272551) q[2];
sx q[2];
rz(-1.9151245) q[2];
sx q[2];
rz(1.2862213) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1133729) q[1];
sx q[1];
rz(-1.7689147) q[1];
sx q[1];
rz(-3.1199725) q[1];
rz(-pi) q[2];
x q[2];
rz(1.008607) q[3];
sx q[3];
rz(-2.6531124) q[3];
sx q[3];
rz(-2.307596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.10716001) q[2];
sx q[2];
rz(-1.1792504) q[2];
sx q[2];
rz(-3.072928) q[2];
rz(-2.5717403) q[3];
sx q[3];
rz(-0.48044258) q[3];
sx q[3];
rz(2.7710052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.079085199) q[0];
sx q[0];
rz(-2.469049) q[0];
sx q[0];
rz(-1.3154718) q[0];
rz(1.8585809) q[1];
sx q[1];
rz(-2.7048769) q[1];
sx q[1];
rz(3.0924996) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21441653) q[0];
sx q[0];
rz(-0.11888725) q[0];
sx q[0];
rz(-2.2772339) q[0];
rz(-pi) q[1];
rz(-2.1911088) q[2];
sx q[2];
rz(-1.1178218) q[2];
sx q[2];
rz(-2.9633303) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9994574) q[1];
sx q[1];
rz(-1.9519567) q[1];
sx q[1];
rz(0.8183523) q[1];
rz(-pi) q[2];
rz(3.0144948) q[3];
sx q[3];
rz(-2.2107901) q[3];
sx q[3];
rz(-0.52332969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7580238) q[2];
sx q[2];
rz(-2.263676) q[2];
sx q[2];
rz(-0.54086584) q[2];
rz(1.0848378) q[3];
sx q[3];
rz(-0.70298755) q[3];
sx q[3];
rz(1.7613523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.908602) q[0];
sx q[0];
rz(-2.1451696) q[0];
sx q[0];
rz(0.10502271) q[0];
rz(0.54166334) q[1];
sx q[1];
rz(-0.88625208) q[1];
sx q[1];
rz(0.35596102) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0070156688) q[0];
sx q[0];
rz(-1.7099713) q[0];
sx q[0];
rz(-0.037875847) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6554285) q[2];
sx q[2];
rz(-0.95900671) q[2];
sx q[2];
rz(-2.391784) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7374266) q[1];
sx q[1];
rz(-1.9475137) q[1];
sx q[1];
rz(-2.1884723) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0186986) q[3];
sx q[3];
rz(-1.8884522) q[3];
sx q[3];
rz(-1.4113219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9291222) q[2];
sx q[2];
rz(-1.4996303) q[2];
sx q[2];
rz(1.8222202) q[2];
rz(-1.3245373) q[3];
sx q[3];
rz(-1.5683441) q[3];
sx q[3];
rz(-0.96778473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20973715) q[0];
sx q[0];
rz(-3.1071438) q[0];
sx q[0];
rz(-1.4631648) q[0];
rz(1.6819008) q[1];
sx q[1];
rz(-1.9821143) q[1];
sx q[1];
rz(-0.34585888) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5992085) q[0];
sx q[0];
rz(-3.0488692) q[0];
sx q[0];
rz(1.0371764) q[0];
rz(-pi) q[1];
rz(1.6357291) q[2];
sx q[2];
rz(-1.254515) q[2];
sx q[2];
rz(0.3102749) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.74259669) q[1];
sx q[1];
rz(-0.95539871) q[1];
sx q[1];
rz(-2.1324498) q[1];
rz(1.1283952) q[3];
sx q[3];
rz(-0.69303362) q[3];
sx q[3];
rz(2.5062163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2807002) q[2];
sx q[2];
rz(-1.8701376) q[2];
sx q[2];
rz(2.8806809) q[2];
rz(-1.4704618) q[3];
sx q[3];
rz(-1.7390395) q[3];
sx q[3];
rz(1.7117333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83229257) q[0];
sx q[0];
rz(-1.1080879) q[0];
sx q[0];
rz(2.9453887) q[0];
rz(0.55083864) q[1];
sx q[1];
rz(-1.6549587) q[1];
sx q[1];
rz(0.5400198) q[1];
rz(0.84340855) q[2];
sx q[2];
rz(-1.7512097) q[2];
sx q[2];
rz(-0.51633121) q[2];
rz(1.2960008) q[3];
sx q[3];
rz(-0.85920371) q[3];
sx q[3];
rz(0.68026713) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
