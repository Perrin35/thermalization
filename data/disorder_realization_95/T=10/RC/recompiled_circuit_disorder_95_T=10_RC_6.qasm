OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.99958217) q[0];
sx q[0];
rz(5.2810623) q[0];
sx q[0];
rz(5.3856344) q[0];
rz(-0.23437962) q[1];
sx q[1];
rz(-0.27581629) q[1];
sx q[1];
rz(2.0770567) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2250741) q[0];
sx q[0];
rz(-1.6406035) q[0];
sx q[0];
rz(2.201447) q[0];
x q[1];
rz(-1.1041553) q[2];
sx q[2];
rz(-0.68744126) q[2];
sx q[2];
rz(1.825037) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7251016) q[1];
sx q[1];
rz(-0.93826586) q[1];
sx q[1];
rz(1.3593258) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3834121) q[3];
sx q[3];
rz(-1.7054134) q[3];
sx q[3];
rz(-1.7077703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.477318) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(-0.73195362) q[2];
rz(0.96015635) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(-1.7094973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9663548) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(-0.91645855) q[0];
rz(-2.6610999) q[1];
sx q[1];
rz(-2.5669211) q[1];
sx q[1];
rz(2.2629471) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71683305) q[0];
sx q[0];
rz(-2.3905907) q[0];
sx q[0];
rz(0.99320937) q[0];
rz(1.8354561) q[2];
sx q[2];
rz(-0.96894962) q[2];
sx q[2];
rz(-0.92483172) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2616927) q[1];
sx q[1];
rz(-2.7163134) q[1];
sx q[1];
rz(-0.17875032) q[1];
rz(-pi) q[2];
rz(-1.8804178) q[3];
sx q[3];
rz(-1.6358346) q[3];
sx q[3];
rz(-2.0464954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2800704) q[2];
sx q[2];
rz(-1.4081988) q[2];
sx q[2];
rz(2.4943165) q[2];
rz(-0.17368008) q[3];
sx q[3];
rz(-2.0208385) q[3];
sx q[3];
rz(2.98996) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6140401) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(2.8955984) q[0];
rz(1.7315158) q[1];
sx q[1];
rz(-1.1743816) q[1];
sx q[1];
rz(-1.1211959) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66378731) q[0];
sx q[0];
rz(-0.87513798) q[0];
sx q[0];
rz(-1.0870766) q[0];
rz(-1.4175225) q[2];
sx q[2];
rz(-2.7445265) q[2];
sx q[2];
rz(-3.0561471) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1998636) q[1];
sx q[1];
rz(-2.6911096) q[1];
sx q[1];
rz(-0.73168879) q[1];
x q[2];
rz(0.21568732) q[3];
sx q[3];
rz(-0.50369278) q[3];
sx q[3];
rz(-1.7220725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.40144172) q[2];
sx q[2];
rz(-1.4462877) q[2];
sx q[2];
rz(2.1901954) q[2];
rz(2.4915063) q[3];
sx q[3];
rz(-1.254046) q[3];
sx q[3];
rz(-0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.6275416) q[0];
sx q[0];
rz(-2.5294332) q[0];
sx q[0];
rz(-0.78805584) q[0];
rz(2.0846941) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(-3.025211) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0375835) q[0];
sx q[0];
rz(-1.5452256) q[0];
sx q[0];
rz(2.5081162) q[0];
rz(-1.0668683) q[2];
sx q[2];
rz(-1.908784) q[2];
sx q[2];
rz(0.68901125) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0105646) q[1];
sx q[1];
rz(-0.66337913) q[1];
sx q[1];
rz(2.4478711) q[1];
rz(-pi) q[2];
rz(-3.0451123) q[3];
sx q[3];
rz(-0.49177846) q[3];
sx q[3];
rz(1.552856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5300166) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(-0.25203618) q[2];
rz(-2.7633372) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(-2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5761121) q[0];
sx q[0];
rz(-1.7385372) q[0];
sx q[0];
rz(-2.6089923) q[0];
rz(-1.4415007) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(-1.1486357) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7589353) q[0];
sx q[0];
rz(-1.4969345) q[0];
sx q[0];
rz(0.9473243) q[0];
x q[1];
rz(-0.45102851) q[2];
sx q[2];
rz(-2.1941059) q[2];
sx q[2];
rz(-0.23652467) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3682813) q[1];
sx q[1];
rz(-1.2174264) q[1];
sx q[1];
rz(1.3694622) q[1];
x q[2];
rz(2.4162021) q[3];
sx q[3];
rz(-1.8728421) q[3];
sx q[3];
rz(-2.6186752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1308412) q[2];
sx q[2];
rz(-0.66117078) q[2];
sx q[2];
rz(2.1095236) q[2];
rz(-0.71470913) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(-0.023795279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59153581) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(2.7456039) q[0];
rz(-1.6962601) q[1];
sx q[1];
rz(-1.6991801) q[1];
sx q[1];
rz(0.2063624) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9128742) q[0];
sx q[0];
rz(-1.4836856) q[0];
sx q[0];
rz(-0.74761439) q[0];
rz(0.80870734) q[2];
sx q[2];
rz(-0.84918298) q[2];
sx q[2];
rz(-1.8546113) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3316162) q[1];
sx q[1];
rz(-0.14252256) q[1];
sx q[1];
rz(-2.317821) q[1];
rz(-pi) q[2];
rz(-2.9456375) q[3];
sx q[3];
rz(-1.2489508) q[3];
sx q[3];
rz(0.67924196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6017194) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(-0.27077857) q[2];
rz(-0.74832908) q[3];
sx q[3];
rz(-1.7691408) q[3];
sx q[3];
rz(1.07871) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2449743) q[0];
sx q[0];
rz(-1.1029607) q[0];
sx q[0];
rz(2.5653429) q[0];
rz(-0.17310625) q[1];
sx q[1];
rz(-0.74179596) q[1];
sx q[1];
rz(1.2111838) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21748397) q[0];
sx q[0];
rz(-2.3586914) q[0];
sx q[0];
rz(-2.3379675) q[0];
rz(-0.25522916) q[2];
sx q[2];
rz(-0.50349871) q[2];
sx q[2];
rz(0.90571678) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2562099) q[1];
sx q[1];
rz(-1.4146283) q[1];
sx q[1];
rz(-0.068710879) q[1];
rz(0.32254036) q[3];
sx q[3];
rz(-2.1493559) q[3];
sx q[3];
rz(1.6998147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8048191) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(0.024519196) q[2];
rz(2.426614) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(-1.582675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1361168) q[0];
sx q[0];
rz(-1.5908717) q[0];
sx q[0];
rz(-1.0205644) q[0];
rz(-0.15469805) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(1.221009) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3295791) q[0];
sx q[0];
rz(-2.3276969) q[0];
sx q[0];
rz(1.194186) q[0];
rz(0.17692716) q[2];
sx q[2];
rz(-1.4104341) q[2];
sx q[2];
rz(-0.71080506) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.80291058) q[1];
sx q[1];
rz(-1.7065587) q[1];
sx q[1];
rz(1.459383) q[1];
rz(-0.33631781) q[3];
sx q[3];
rz(-0.85018966) q[3];
sx q[3];
rz(-1.1433126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.28875479) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(0.70518804) q[2];
rz(-2.5382036) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(-1.6220629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.110638) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(3.0045793) q[0];
rz(2.5367472) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(-0.12577122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98709479) q[0];
sx q[0];
rz(-0.25941601) q[0];
sx q[0];
rz(0.82018606) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.618082) q[2];
sx q[2];
rz(-2.7611809) q[2];
sx q[2];
rz(-1.7720122) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3198493) q[1];
sx q[1];
rz(-0.84801596) q[1];
sx q[1];
rz(-2.7601932) q[1];
x q[2];
rz(-0.59154193) q[3];
sx q[3];
rz(-0.35174832) q[3];
sx q[3];
rz(1.2951617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.003309) q[2];
sx q[2];
rz(-0.97390276) q[2];
sx q[2];
rz(-0.70927817) q[2];
rz(-0.62018958) q[3];
sx q[3];
rz(-0.87738335) q[3];
sx q[3];
rz(2.9848849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9897292) q[0];
sx q[0];
rz(-2.3487838) q[0];
sx q[0];
rz(1.2003157) q[0];
rz(2.8740846) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(-1.1402003) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4338019) q[0];
sx q[0];
rz(-1.8330488) q[0];
sx q[0];
rz(-3.0387525) q[0];
rz(-2.8814949) q[2];
sx q[2];
rz(-2.0195228) q[2];
sx q[2];
rz(1.5425494) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4789341) q[1];
sx q[1];
rz(-1.3824029) q[1];
sx q[1];
rz(-1.4807832) q[1];
x q[2];
rz(-0.31302281) q[3];
sx q[3];
rz(-0.81837294) q[3];
sx q[3];
rz(2.5767874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7662979) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(0.89861384) q[2];
rz(-3.1344154) q[3];
sx q[3];
rz(-2.3982748) q[3];
sx q[3];
rz(1.7858508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2469149) q[0];
sx q[0];
rz(-1.4151731) q[0];
sx q[0];
rz(-1.7807501) q[0];
rz(-2.6208411) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(-1.0976787) q[2];
sx q[2];
rz(-1.1259148) q[2];
sx q[2];
rz(-1.7150707) q[2];
rz(-2.5614212) q[3];
sx q[3];
rz(-0.67065722) q[3];
sx q[3];
rz(-1.8967659) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
