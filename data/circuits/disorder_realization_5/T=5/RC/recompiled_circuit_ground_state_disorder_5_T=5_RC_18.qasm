OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4910645) q[0];
sx q[0];
rz(-2.129038) q[0];
sx q[0];
rz(-0.92232409) q[0];
rz(-0.84696472) q[1];
sx q[1];
rz(-1.6672517) q[1];
sx q[1];
rz(0.21811952) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51061741) q[0];
sx q[0];
rz(-1.4790863) q[0];
sx q[0];
rz(2.6789078) q[0];
rz(-pi) q[1];
rz(-1.5297024) q[2];
sx q[2];
rz(-1.1930704) q[2];
sx q[2];
rz(-0.95313493) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3936945) q[1];
sx q[1];
rz(-1.6130578) q[1];
sx q[1];
rz(1.7903922) q[1];
x q[2];
rz(-1.140959) q[3];
sx q[3];
rz(-1.6419193) q[3];
sx q[3];
rz(1.5971368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.073079022) q[2];
sx q[2];
rz(-1.928669) q[2];
sx q[2];
rz(2.6050513) q[2];
rz(-1.0088629) q[3];
sx q[3];
rz(-0.77102414) q[3];
sx q[3];
rz(2.3276276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5040078) q[0];
sx q[0];
rz(-1.3115839) q[0];
sx q[0];
rz(3.1245533) q[0];
rz(-1.0631961) q[1];
sx q[1];
rz(-0.76779643) q[1];
sx q[1];
rz(-2.1868736) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9794036) q[0];
sx q[0];
rz(-1.5155751) q[0];
sx q[0];
rz(0.22947854) q[0];
x q[1];
rz(1.0771712) q[2];
sx q[2];
rz(-0.33702484) q[2];
sx q[2];
rz(0.21060196) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0481457) q[1];
sx q[1];
rz(-2.8110782) q[1];
sx q[1];
rz(-2.3814209) q[1];
rz(-pi) q[2];
rz(-1.6240261) q[3];
sx q[3];
rz(-1.551515) q[3];
sx q[3];
rz(-0.35010168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9118328) q[2];
sx q[2];
rz(-2.2288897) q[2];
sx q[2];
rz(-1.1141874) q[2];
rz(-2.0071425) q[3];
sx q[3];
rz(-2.9454234) q[3];
sx q[3];
rz(0.1964868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.590362) q[0];
sx q[0];
rz(-0.95142618) q[0];
sx q[0];
rz(2.9587342) q[0];
rz(-0.081347801) q[1];
sx q[1];
rz(-2.3902939) q[1];
sx q[1];
rz(2.523211) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30899042) q[0];
sx q[0];
rz(-1.6438577) q[0];
sx q[0];
rz(-2.354524) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19962991) q[2];
sx q[2];
rz(-0.77130328) q[2];
sx q[2];
rz(-0.16801258) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.394819) q[1];
sx q[1];
rz(-1.2319984) q[1];
sx q[1];
rz(2.9940786) q[1];
rz(2.5181055) q[3];
sx q[3];
rz(-2.1716431) q[3];
sx q[3];
rz(-2.2740325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1387834) q[2];
sx q[2];
rz(-1.8178136) q[2];
sx q[2];
rz(0.36661026) q[2];
rz(-1.2159411) q[3];
sx q[3];
rz(-1.9462908) q[3];
sx q[3];
rz(-0.6634357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4902896) q[0];
sx q[0];
rz(-1.2325352) q[0];
sx q[0];
rz(-2.41462) q[0];
rz(-0.90826774) q[1];
sx q[1];
rz(-0.5779225) q[1];
sx q[1];
rz(1.04331) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7730966) q[0];
sx q[0];
rz(-0.90963042) q[0];
sx q[0];
rz(-1.5776683) q[0];
rz(-pi) q[1];
rz(1.415148) q[2];
sx q[2];
rz(-0.37964941) q[2];
sx q[2];
rz(1.8774892) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9227495) q[1];
sx q[1];
rz(-2.6237539) q[1];
sx q[1];
rz(-1.0475558) q[1];
rz(-pi) q[2];
rz(-2.7611012) q[3];
sx q[3];
rz(-1.0416608) q[3];
sx q[3];
rz(2.7545415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4116481) q[2];
sx q[2];
rz(-0.32117716) q[2];
sx q[2];
rz(1.8357065) q[2];
rz(-3.0236687) q[3];
sx q[3];
rz(-2.6730461) q[3];
sx q[3];
rz(-2.0809295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0333198) q[0];
sx q[0];
rz(-0.98618996) q[0];
sx q[0];
rz(1.6075851) q[0];
rz(2.8458505) q[1];
sx q[1];
rz(-1.60138) q[1];
sx q[1];
rz(1.6827513) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0223607) q[0];
sx q[0];
rz(-1.5007449) q[0];
sx q[0];
rz(2.3168867) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2590911) q[2];
sx q[2];
rz(-2.6771328) q[2];
sx q[2];
rz(-1.8732173) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0572059) q[1];
sx q[1];
rz(-0.73250341) q[1];
sx q[1];
rz(-1.8671579) q[1];
x q[2];
rz(-1.9146054) q[3];
sx q[3];
rz(-0.60887486) q[3];
sx q[3];
rz(-0.46907779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7532588) q[2];
sx q[2];
rz(-0.39187852) q[2];
sx q[2];
rz(-0.64588109) q[2];
rz(-1.9258457) q[3];
sx q[3];
rz(-1.4589717) q[3];
sx q[3];
rz(-1.7997883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011768613) q[0];
sx q[0];
rz(-2.3079066) q[0];
sx q[0];
rz(-2.5339793) q[0];
rz(-1.1786849) q[1];
sx q[1];
rz(-1.8557529) q[1];
sx q[1];
rz(-2.7778621) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0775438) q[0];
sx q[0];
rz(-2.0850967) q[0];
sx q[0];
rz(-0.27405996) q[0];
rz(-pi) q[1];
rz(-1.7044439) q[2];
sx q[2];
rz(-2.8060348) q[2];
sx q[2];
rz(0.37119532) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2567098) q[1];
sx q[1];
rz(-2.3824168) q[1];
sx q[1];
rz(-1.3961387) q[1];
rz(-pi) q[2];
rz(-2.167114) q[3];
sx q[3];
rz(-1.5913561) q[3];
sx q[3];
rz(-2.1316776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.58934775) q[2];
sx q[2];
rz(-3.037368) q[2];
sx q[2];
rz(-1.1927401) q[2];
rz(-0.73181152) q[3];
sx q[3];
rz(-1.4251499) q[3];
sx q[3];
rz(1.4881136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8530387) q[0];
sx q[0];
rz(-2.9712501) q[0];
sx q[0];
rz(-1.6188251) q[0];
rz(1.4672) q[1];
sx q[1];
rz(-1.3102945) q[1];
sx q[1];
rz(0.67749643) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8554133) q[0];
sx q[0];
rz(-0.87550801) q[0];
sx q[0];
rz(3.0762818) q[0];
rz(-0.23394312) q[2];
sx q[2];
rz(-1.3877227) q[2];
sx q[2];
rz(1.853136) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.575857) q[1];
sx q[1];
rz(-1.1850372) q[1];
sx q[1];
rz(1.1123453) q[1];
rz(-2.4091408) q[3];
sx q[3];
rz(-1.4518514) q[3];
sx q[3];
rz(0.67057395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39685321) q[2];
sx q[2];
rz(-0.91529673) q[2];
sx q[2];
rz(1.3057365) q[2];
rz(2.8030677) q[3];
sx q[3];
rz(-1.1208231) q[3];
sx q[3];
rz(-2.4566076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6029538) q[0];
sx q[0];
rz(-0.21553497) q[0];
sx q[0];
rz(0.18390528) q[0];
rz(1.6869847) q[1];
sx q[1];
rz(-0.55211663) q[1];
sx q[1];
rz(-0.49585453) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5372779) q[0];
sx q[0];
rz(-2.0576982) q[0];
sx q[0];
rz(-0.30211289) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3073436) q[2];
sx q[2];
rz(-0.96820346) q[2];
sx q[2];
rz(-1.955223) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.293893) q[1];
sx q[1];
rz(-1.9260848) q[1];
sx q[1];
rz(-1.9778663) q[1];
rz(-pi) q[2];
rz(2.1231648) q[3];
sx q[3];
rz(-1.4694858) q[3];
sx q[3];
rz(0.61063572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0869861) q[2];
sx q[2];
rz(-2.2973674) q[2];
sx q[2];
rz(1.9795214) q[2];
rz(-1.453513) q[3];
sx q[3];
rz(-0.96265692) q[3];
sx q[3];
rz(-1.2801722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7206955) q[0];
sx q[0];
rz(-1.9891885) q[0];
sx q[0];
rz(-2.5757117) q[0];
rz(0.3903009) q[1];
sx q[1];
rz(-1.5696328) q[1];
sx q[1];
rz(1.3077259) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54057656) q[0];
sx q[0];
rz(-1.8800288) q[0];
sx q[0];
rz(1.0544257) q[0];
rz(-pi) q[1];
rz(2.0686223) q[2];
sx q[2];
rz(-2.7171869) q[2];
sx q[2];
rz(2.5967732) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81586397) q[1];
sx q[1];
rz(-0.96885175) q[1];
sx q[1];
rz(-1.8900327) q[1];
x q[2];
rz(2.2447137) q[3];
sx q[3];
rz(-0.19052902) q[3];
sx q[3];
rz(0.49303699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.86964837) q[2];
sx q[2];
rz(-2.7248236) q[2];
sx q[2];
rz(-2.1040253) q[2];
rz(-1.3197673) q[3];
sx q[3];
rz(-0.82873738) q[3];
sx q[3];
rz(2.493609) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8808402) q[0];
sx q[0];
rz(-0.97761959) q[0];
sx q[0];
rz(-2.4993437) q[0];
rz(1.3308659) q[1];
sx q[1];
rz(-0.86527491) q[1];
sx q[1];
rz(-0.16255249) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2519788) q[0];
sx q[0];
rz(-1.9601213) q[0];
sx q[0];
rz(-2.2963803) q[0];
rz(-2.3847975) q[2];
sx q[2];
rz(-0.50163236) q[2];
sx q[2];
rz(-1.0197786) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6813415) q[1];
sx q[1];
rz(-0.5149018) q[1];
sx q[1];
rz(-1.9327362) q[1];
rz(-pi) q[2];
rz(1.4914042) q[3];
sx q[3];
rz(-2.1526436) q[3];
sx q[3];
rz(1.013085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.51005298) q[2];
sx q[2];
rz(-1.2756462) q[2];
sx q[2];
rz(-2.1007382) q[2];
rz(-0.96796525) q[3];
sx q[3];
rz(-1.0699882) q[3];
sx q[3];
rz(0.66106838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80617245) q[0];
sx q[0];
rz(-1.5652884) q[0];
sx q[0];
rz(-0.096927222) q[0];
rz(-2.9087635) q[1];
sx q[1];
rz(-1.0284582) q[1];
sx q[1];
rz(-3.1340541) q[1];
rz(2.0816879) q[2];
sx q[2];
rz(-2.2507406) q[2];
sx q[2];
rz(3.0377664) q[2];
rz(2.4416853) q[3];
sx q[3];
rz(-1.2716765) q[3];
sx q[3];
rz(-1.8636462) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
