OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7656443) q[0];
sx q[0];
rz(-0.41416895) q[0];
sx q[0];
rz(2.3400657) q[0];
rz(-0.0030567788) q[1];
sx q[1];
rz(-0.86441511) q[1];
sx q[1];
rz(-0.094223082) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021359062) q[0];
sx q[0];
rz(-2.2968043) q[0];
sx q[0];
rz(-2.0840706) q[0];
x q[1];
rz(-1.8434486) q[2];
sx q[2];
rz(-1.6003719) q[2];
sx q[2];
rz(0.77405133) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7677418) q[1];
sx q[1];
rz(-2.5818425) q[1];
sx q[1];
rz(2.1327726) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5729957) q[3];
sx q[3];
rz(-1.723395) q[3];
sx q[3];
rz(-0.27657498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6039383) q[2];
sx q[2];
rz(-2.5014169) q[2];
sx q[2];
rz(1.5113277) q[2];
rz(2.9553735) q[3];
sx q[3];
rz(-2.7112466) q[3];
sx q[3];
rz(-0.65096861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6179825) q[0];
sx q[0];
rz(-1.9600493) q[0];
sx q[0];
rz(0.13667662) q[0];
rz(-3.0467721) q[1];
sx q[1];
rz(-2.6780728) q[1];
sx q[1];
rz(-0.76900855) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.626132) q[0];
sx q[0];
rz(-2.7410586) q[0];
sx q[0];
rz(-0.35450165) q[0];
rz(-1.9273754) q[2];
sx q[2];
rz(-0.90889034) q[2];
sx q[2];
rz(1.8476768) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5040759) q[1];
sx q[1];
rz(-0.83607644) q[1];
sx q[1];
rz(-1.3512263) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7058825) q[3];
sx q[3];
rz(-1.6052433) q[3];
sx q[3];
rz(-2.5273539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7424221) q[2];
sx q[2];
rz(-2.6593282) q[2];
sx q[2];
rz(1.4313618) q[2];
rz(0.4906022) q[3];
sx q[3];
rz(-0.49121818) q[3];
sx q[3];
rz(2.365999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8785777) q[0];
sx q[0];
rz(-2.7624625) q[0];
sx q[0];
rz(-1.0947134) q[0];
rz(1.3021944) q[1];
sx q[1];
rz(-2.1523988) q[1];
sx q[1];
rz(3.1375695) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99056315) q[0];
sx q[0];
rz(-0.73712611) q[0];
sx q[0];
rz(1.766979) q[0];
rz(-pi) q[1];
rz(-2.4039359) q[2];
sx q[2];
rz(-1.0152384) q[2];
sx q[2];
rz(0.45762024) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.04847542) q[1];
sx q[1];
rz(-0.25893908) q[1];
sx q[1];
rz(3.0228651) q[1];
rz(-1.324292) q[3];
sx q[3];
rz(-0.87216264) q[3];
sx q[3];
rz(2.7072631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5502988) q[2];
sx q[2];
rz(-2.55547) q[2];
sx q[2];
rz(2.6934521) q[2];
rz(2.6552933) q[3];
sx q[3];
rz(-0.86217642) q[3];
sx q[3];
rz(1.1264616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9124741) q[0];
sx q[0];
rz(-1.7121226) q[0];
sx q[0];
rz(0.66977704) q[0];
rz(-1.229333) q[1];
sx q[1];
rz(-0.45893097) q[1];
sx q[1];
rz(-2.5965447) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3236967) q[0];
sx q[0];
rz(-2.4113089) q[0];
sx q[0];
rz(-0.72636168) q[0];
rz(2.2322234) q[2];
sx q[2];
rz(-2.0880359) q[2];
sx q[2];
rz(1.478372) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4008444) q[1];
sx q[1];
rz(-1.6312338) q[1];
sx q[1];
rz(1.7857741) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7983851) q[3];
sx q[3];
rz(-1.9163864) q[3];
sx q[3];
rz(-1.8963199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5547319) q[2];
sx q[2];
rz(-1.6394337) q[2];
sx q[2];
rz(-0.16793212) q[2];
rz(-1.2083017) q[3];
sx q[3];
rz(-0.30253634) q[3];
sx q[3];
rz(-1.6944073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55030441) q[0];
sx q[0];
rz(-1.8809603) q[0];
sx q[0];
rz(3.1377129) q[0];
rz(-0.30715352) q[1];
sx q[1];
rz(-2.3852564) q[1];
sx q[1];
rz(-0.53877962) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2821747) q[0];
sx q[0];
rz(-0.27056405) q[0];
sx q[0];
rz(2.5859588) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4319375) q[2];
sx q[2];
rz(-1.1449688) q[2];
sx q[2];
rz(2.6151997) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6314124) q[1];
sx q[1];
rz(-2.4459527) q[1];
sx q[1];
rz(-2.5744777) q[1];
rz(-2.0224376) q[3];
sx q[3];
rz(-2.2502101) q[3];
sx q[3];
rz(1.7070626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6946081) q[2];
sx q[2];
rz(-2.0422523) q[2];
sx q[2];
rz(1.851932) q[2];
rz(1.6192216) q[3];
sx q[3];
rz(-0.74519849) q[3];
sx q[3];
rz(1.4815909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8412142) q[0];
sx q[0];
rz(-2.2024246) q[0];
sx q[0];
rz(-0.10264957) q[0];
rz(2.4094021) q[1];
sx q[1];
rz(-2.3468572) q[1];
sx q[1];
rz(0.57714677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2271778) q[0];
sx q[0];
rz(-0.028870048) q[0];
sx q[0];
rz(-1.6976852) q[0];
rz(-2.5546165) q[2];
sx q[2];
rz(-1.7086626) q[2];
sx q[2];
rz(-1.1816813) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.53458339) q[1];
sx q[1];
rz(-1.4614236) q[1];
sx q[1];
rz(-1.4916199) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8982557) q[3];
sx q[3];
rz(-0.84210448) q[3];
sx q[3];
rz(-3.1407243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.010043667) q[2];
sx q[2];
rz(-2.4250344) q[2];
sx q[2];
rz(-1.6183759) q[2];
rz(3.0736382) q[3];
sx q[3];
rz(-0.32502919) q[3];
sx q[3];
rz(-2.3006191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0129358) q[0];
sx q[0];
rz(-2.1067297) q[0];
sx q[0];
rz(0.77142429) q[0];
rz(-2.6029288) q[1];
sx q[1];
rz(-1.795105) q[1];
sx q[1];
rz(-1.0661941) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5911884) q[0];
sx q[0];
rz(-0.26445358) q[0];
sx q[0];
rz(1.868684) q[0];
rz(-pi) q[1];
rz(0.23287878) q[2];
sx q[2];
rz(-1.0693197) q[2];
sx q[2];
rz(-1.6362783) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.33246189) q[1];
sx q[1];
rz(-1.5576502) q[1];
sx q[1];
rz(2.9244208) q[1];
x q[2];
rz(-0.33447075) q[3];
sx q[3];
rz(-1.3190509) q[3];
sx q[3];
rz(-1.9289964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.44925877) q[2];
sx q[2];
rz(-0.80654311) q[2];
sx q[2];
rz(-2.6991357) q[2];
rz(2.9509406) q[3];
sx q[3];
rz(-0.72398829) q[3];
sx q[3];
rz(0.75907069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-2.8964748) q[0];
sx q[0];
rz(-2.9560095) q[0];
sx q[0];
rz(1.3099439) q[0];
rz(-2.7438121) q[1];
sx q[1];
rz(-2.3186627) q[1];
sx q[1];
rz(-1.3936874) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8373972) q[0];
sx q[0];
rz(-1.0365573) q[0];
sx q[0];
rz(2.6287352) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23598052) q[2];
sx q[2];
rz(-1.2771291) q[2];
sx q[2];
rz(-2.8666509) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9164394) q[1];
sx q[1];
rz(-1.272456) q[1];
sx q[1];
rz(-2.6544996) q[1];
x q[2];
rz(-1.769479) q[3];
sx q[3];
rz(-1.6749089) q[3];
sx q[3];
rz(2.7840691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1132249) q[2];
sx q[2];
rz(-2.4312879) q[2];
sx q[2];
rz(-0.073471546) q[2];
rz(2.4469589) q[3];
sx q[3];
rz(-1.8690542) q[3];
sx q[3];
rz(-1.0276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46591127) q[0];
sx q[0];
rz(-0.2008734) q[0];
sx q[0];
rz(0.95941108) q[0];
rz(-0.77964669) q[1];
sx q[1];
rz(-2.1387073) q[1];
sx q[1];
rz(2.8693105) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.773943) q[0];
sx q[0];
rz(-3.127357) q[0];
sx q[0];
rz(-1.9444501) q[0];
rz(-pi) q[1];
x q[1];
rz(2.831976) q[2];
sx q[2];
rz(-0.88954347) q[2];
sx q[2];
rz(2.3818784) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4466557) q[1];
sx q[1];
rz(-0.65498252) q[1];
sx q[1];
rz(0.86465624) q[1];
rz(-pi) q[2];
rz(0.7471146) q[3];
sx q[3];
rz(-1.1076339) q[3];
sx q[3];
rz(0.85827352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4043364) q[2];
sx q[2];
rz(-1.6360838) q[2];
sx q[2];
rz(-0.20757248) q[2];
rz(0.10457822) q[3];
sx q[3];
rz(-0.50555491) q[3];
sx q[3];
rz(-1.6149717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614457) q[0];
sx q[0];
rz(-1.9141645) q[0];
sx q[0];
rz(-2.0848059) q[0];
rz(2.5959004) q[1];
sx q[1];
rz(-0.97815198) q[1];
sx q[1];
rz(-0.32870764) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0595039) q[0];
sx q[0];
rz(-0.57841136) q[0];
sx q[0];
rz(-2.2057981) q[0];
rz(-pi) q[1];
rz(1.516009) q[2];
sx q[2];
rz(-1.6431) q[2];
sx q[2];
rz(-0.26501467) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4663965) q[1];
sx q[1];
rz(-1.5874784) q[1];
sx q[1];
rz(1.9128602) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8449109) q[3];
sx q[3];
rz(-2.3233827) q[3];
sx q[3];
rz(-0.75358281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1591961) q[2];
sx q[2];
rz(-3.0089162) q[2];
sx q[2];
rz(-1.1245493) q[2];
rz(3.0647965) q[3];
sx q[3];
rz(-0.43282893) q[3];
sx q[3];
rz(2.2176149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9818253) q[0];
sx q[0];
rz(-1.8851017) q[0];
sx q[0];
rz(-1.5633352) q[0];
rz(0.81623296) q[1];
sx q[1];
rz(-1.4030133) q[1];
sx q[1];
rz(-1.0907008) q[1];
rz(-2.6098567) q[2];
sx q[2];
rz(-2.3219863) q[2];
sx q[2];
rz(1.815997) q[2];
rz(3.1327457) q[3];
sx q[3];
rz(-2.736241) q[3];
sx q[3];
rz(-0.74373087) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
