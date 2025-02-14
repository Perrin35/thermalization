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
rz(-0.24902046) q[0];
sx q[0];
rz(-0.95872107) q[0];
sx q[0];
rz(-0.22554654) q[0];
rz(-1.072999) q[1];
sx q[1];
rz(3.5993242) q[1];
sx q[1];
rz(9.2158894) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6207558) q[0];
sx q[0];
rz(-2.0055456) q[0];
sx q[0];
rz(-1.2723544) q[0];
rz(-pi) q[1];
rz(1.289008) q[2];
sx q[2];
rz(-0.95880858) q[2];
sx q[2];
rz(-2.2300697) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35395216) q[1];
sx q[1];
rz(-0.6122511) q[1];
sx q[1];
rz(2.5933215) q[1];
rz(-pi) q[2];
rz(2.3129102) q[3];
sx q[3];
rz(-0.4503612) q[3];
sx q[3];
rz(-0.27638232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0245584) q[2];
sx q[2];
rz(-0.25091761) q[2];
sx q[2];
rz(-0.92570242) q[2];
rz(0.77445817) q[3];
sx q[3];
rz(-1.9175994) q[3];
sx q[3];
rz(-1.8584049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16933146) q[0];
sx q[0];
rz(-0.65003482) q[0];
sx q[0];
rz(1.6768804) q[0];
rz(0.88893923) q[1];
sx q[1];
rz(-2.2496532) q[1];
sx q[1];
rz(-1.3533786) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66799712) q[0];
sx q[0];
rz(-1.8596853) q[0];
sx q[0];
rz(2.5210565) q[0];
rz(-pi) q[1];
rz(0.057463138) q[2];
sx q[2];
rz(-2.7749535) q[2];
sx q[2];
rz(-2.5861248) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6135237) q[1];
sx q[1];
rz(-2.1557243) q[1];
sx q[1];
rz(1.3207573) q[1];
x q[2];
rz(-2.3634745) q[3];
sx q[3];
rz(-1.3740842) q[3];
sx q[3];
rz(-0.16315483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1464403) q[2];
sx q[2];
rz(-2.1285987) q[2];
sx q[2];
rz(2.1688482) q[2];
rz(1.815833) q[3];
sx q[3];
rz(-1.1498068) q[3];
sx q[3];
rz(0.61746922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4165118) q[0];
sx q[0];
rz(-1.4840115) q[0];
sx q[0];
rz(-0.02956477) q[0];
rz(2.120453) q[1];
sx q[1];
rz(-2.857326) q[1];
sx q[1];
rz(2.8254642) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2293623) q[0];
sx q[0];
rz(-1.5490313) q[0];
sx q[0];
rz(-0.020177186) q[0];
x q[1];
rz(-0.093806819) q[2];
sx q[2];
rz(-1.7711763) q[2];
sx q[2];
rz(-2.0721638) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30313213) q[1];
sx q[1];
rz(-0.37835281) q[1];
sx q[1];
rz(-2.0580473) q[1];
rz(-pi) q[2];
rz(-0.088313266) q[3];
sx q[3];
rz(-2.4188015) q[3];
sx q[3];
rz(-0.73107728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.11008392) q[2];
sx q[2];
rz(-1.8434593) q[2];
sx q[2];
rz(0.52424866) q[2];
rz(2.5414069) q[3];
sx q[3];
rz(-0.46349183) q[3];
sx q[3];
rz(2.0704827) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6780739) q[0];
sx q[0];
rz(-1.9215895) q[0];
sx q[0];
rz(-2.779261) q[0];
rz(-2.433297) q[1];
sx q[1];
rz(-0.34168044) q[1];
sx q[1];
rz(-1.1452311) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4701242) q[0];
sx q[0];
rz(-3.0076179) q[0];
sx q[0];
rz(-2.2132316) q[0];
rz(-pi) q[1];
rz(-0.48787222) q[2];
sx q[2];
rz(-0.93624253) q[2];
sx q[2];
rz(1.799859) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.501579) q[1];
sx q[1];
rz(-2.3928436) q[1];
sx q[1];
rz(-1.1925936) q[1];
rz(-pi) q[2];
rz(1.5208552) q[3];
sx q[3];
rz(-1.4068713) q[3];
sx q[3];
rz(2.0633179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.5273253) q[2];
sx q[2];
rz(-2.3136487) q[2];
sx q[2];
rz(3.1329727) q[2];
rz(-2.4873867) q[3];
sx q[3];
rz(-0.71627408) q[3];
sx q[3];
rz(-0.27230984) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9006186) q[0];
sx q[0];
rz(-2.8130377) q[0];
sx q[0];
rz(0.7884489) q[0];
rz(-1.0139326) q[1];
sx q[1];
rz(-1.8221816) q[1];
sx q[1];
rz(-3.0533155) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81472003) q[0];
sx q[0];
rz(-0.86866394) q[0];
sx q[0];
rz(2.4725295) q[0];
rz(2.9806251) q[2];
sx q[2];
rz(-2.5720398) q[2];
sx q[2];
rz(0.19191027) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8535117) q[1];
sx q[1];
rz(-1.8520903) q[1];
sx q[1];
rz(2.0916677) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0658947) q[3];
sx q[3];
rz(-1.5468239) q[3];
sx q[3];
rz(0.8296465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2225515) q[2];
sx q[2];
rz(-1.7307948) q[2];
sx q[2];
rz(2.6338573) q[2];
rz(3.0159085) q[3];
sx q[3];
rz(-1.1832184) q[3];
sx q[3];
rz(-2.3465033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8102201) q[0];
sx q[0];
rz(-0.75646821) q[0];
sx q[0];
rz(3.0561225) q[0];
rz(0.62689176) q[1];
sx q[1];
rz(-1.1566999) q[1];
sx q[1];
rz(-1.8524106) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7580494) q[0];
sx q[0];
rz(-1.6746759) q[0];
sx q[0];
rz(-1.3707177) q[0];
rz(-1.4165381) q[2];
sx q[2];
rz(-0.83310328) q[2];
sx q[2];
rz(1.2546033) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.35297817) q[1];
sx q[1];
rz(-1.7403688) q[1];
sx q[1];
rz(-1.8119078) q[1];
rz(-pi) q[2];
rz(-3.1269642) q[3];
sx q[3];
rz(-1.808015) q[3];
sx q[3];
rz(-0.62495172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36279303) q[2];
sx q[2];
rz(-0.94393602) q[2];
sx q[2];
rz(1.2550521) q[2];
rz(3.0826027) q[3];
sx q[3];
rz(-1.9374282) q[3];
sx q[3];
rz(-2.1571933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7207709) q[0];
sx q[0];
rz(-1.868792) q[0];
sx q[0];
rz(2.3765748) q[0];
rz(-1.7401241) q[1];
sx q[1];
rz(-2.2547289) q[1];
sx q[1];
rz(-2.9041362) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9119806) q[0];
sx q[0];
rz(-1.53541) q[0];
sx q[0];
rz(-2.1982212) q[0];
rz(0.89037322) q[2];
sx q[2];
rz(-1.6294663) q[2];
sx q[2];
rz(0.49505297) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4478893) q[1];
sx q[1];
rz(-2.8475326) q[1];
sx q[1];
rz(2.990688) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2672784) q[3];
sx q[3];
rz(-1.5477936) q[3];
sx q[3];
rz(1.5222331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6852297) q[2];
sx q[2];
rz(-2.3380029) q[2];
sx q[2];
rz(2.2881499) q[2];
rz(1.8026836) q[3];
sx q[3];
rz(-1.5693376) q[3];
sx q[3];
rz(2.9760823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45561871) q[0];
sx q[0];
rz(-1.2116665) q[0];
sx q[0];
rz(-2.9562505) q[0];
rz(2.7475884) q[1];
sx q[1];
rz(-1.3961671) q[1];
sx q[1];
rz(-2.0448763) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1117437) q[0];
sx q[0];
rz(-1.1858479) q[0];
sx q[0];
rz(-2.0864331) q[0];
x q[1];
rz(2.2359879) q[2];
sx q[2];
rz(-1.6597372) q[2];
sx q[2];
rz(0.62155089) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.85307676) q[1];
sx q[1];
rz(-2.4429053) q[1];
sx q[1];
rz(-0.53650155) q[1];
x q[2];
rz(-2.5470252) q[3];
sx q[3];
rz(-1.2487186) q[3];
sx q[3];
rz(-0.26628029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3700221) q[2];
sx q[2];
rz(-3.0574419) q[2];
sx q[2];
rz(2.8328075) q[2];
rz(-1.8874946) q[3];
sx q[3];
rz(-1.2718688) q[3];
sx q[3];
rz(-1.1312283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9517188) q[0];
sx q[0];
rz(-1.254344) q[0];
sx q[0];
rz(0.21981123) q[0];
rz(1.1687357) q[1];
sx q[1];
rz(-1.3558148) q[1];
sx q[1];
rz(-1.2723602) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2975966) q[0];
sx q[0];
rz(-1.5463313) q[0];
sx q[0];
rz(-2.1750227) q[0];
rz(2.6307893) q[2];
sx q[2];
rz(-1.496721) q[2];
sx q[2];
rz(1.6965716) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0366401) q[1];
sx q[1];
rz(-0.4401686) q[1];
sx q[1];
rz(2.7326581) q[1];
rz(1.06274) q[3];
sx q[3];
rz(-1.1241759) q[3];
sx q[3];
rz(1.6376709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9610338) q[2];
sx q[2];
rz(-0.84332931) q[2];
sx q[2];
rz(-2.6832704) q[2];
rz(-0.35442963) q[3];
sx q[3];
rz(-1.7731881) q[3];
sx q[3];
rz(1.9663158) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71427041) q[0];
sx q[0];
rz(-1.9388119) q[0];
sx q[0];
rz(0.74139968) q[0];
rz(1.7085913) q[1];
sx q[1];
rz(-2.7129136) q[1];
sx q[1];
rz(-3.0523849) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7767255) q[0];
sx q[0];
rz(-1.7517543) q[0];
sx q[0];
rz(1.0904161) q[0];
x q[1];
rz(-1.9544425) q[2];
sx q[2];
rz(-2.128278) q[2];
sx q[2];
rz(-1.5889744) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7877823) q[1];
sx q[1];
rz(-1.612376) q[1];
sx q[1];
rz(-1.9099243) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44762917) q[3];
sx q[3];
rz(-2.3530686) q[3];
sx q[3];
rz(-0.41390362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0006813) q[2];
sx q[2];
rz(-0.050364308) q[2];
sx q[2];
rz(-0.64765206) q[2];
rz(-2.7378313) q[3];
sx q[3];
rz(-1.8881366) q[3];
sx q[3];
rz(3.0103179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1972926) q[0];
sx q[0];
rz(-0.97351749) q[0];
sx q[0];
rz(1.0607006) q[0];
rz(2.3464959) q[1];
sx q[1];
rz(-2.2207694) q[1];
sx q[1];
rz(2.3650852) q[1];
rz(0.9153824) q[2];
sx q[2];
rz(-1.6554828) q[2];
sx q[2];
rz(0.2632904) q[2];
rz(-2.285801) q[3];
sx q[3];
rz(-3.0223911) q[3];
sx q[3];
rz(0.08798616) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
