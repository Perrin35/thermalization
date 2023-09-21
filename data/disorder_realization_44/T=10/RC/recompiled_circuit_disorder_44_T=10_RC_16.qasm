OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(6.7232806) q[0];
sx q[0];
rz(6.4203782) q[0];
rz(1.4057012) q[1];
sx q[1];
rz(4.5448137) q[1];
sx q[1];
rz(9.9546976) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7007028) q[0];
sx q[0];
rz(-0.39632495) q[0];
sx q[0];
rz(1.2896145) q[0];
rz(-pi) q[1];
x q[1];
rz(2.797384) q[2];
sx q[2];
rz(-1.1905626) q[2];
sx q[2];
rz(-2.3772079) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6119541) q[1];
sx q[1];
rz(-2.0679592) q[1];
sx q[1];
rz(-2.0938718) q[1];
rz(-pi) q[2];
rz(-0.75833851) q[3];
sx q[3];
rz(-1.3876649) q[3];
sx q[3];
rz(-1.6299712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4522176) q[2];
sx q[2];
rz(-1.8415425) q[2];
sx q[2];
rz(-0.33660647) q[2];
rz(-1.6254788) q[3];
sx q[3];
rz(-0.55364004) q[3];
sx q[3];
rz(-1.5256933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7933554) q[0];
sx q[0];
rz(-2.0331148) q[0];
sx q[0];
rz(-0.021214699) q[0];
rz(-1.1938098) q[1];
sx q[1];
rz(-1.0394916) q[1];
sx q[1];
rz(2.3056727) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6244038) q[0];
sx q[0];
rz(-1.4319341) q[0];
sx q[0];
rz(0.018272321) q[0];
rz(-pi) q[1];
rz(-1.0216653) q[2];
sx q[2];
rz(-1.9252535) q[2];
sx q[2];
rz(-2.3159388) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4513431) q[1];
sx q[1];
rz(-0.93938821) q[1];
sx q[1];
rz(0.33957014) q[1];
rz(-pi) q[2];
rz(-0.73987506) q[3];
sx q[3];
rz(-1.2736819) q[3];
sx q[3];
rz(-1.5418996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8643643) q[2];
sx q[2];
rz(-1.9936864) q[2];
sx q[2];
rz(-1.345984) q[2];
rz(2.7820382) q[3];
sx q[3];
rz(-2.1988726) q[3];
sx q[3];
rz(-2.6446222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42831746) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(-1.0536449) q[0];
rz(-1.2288278) q[1];
sx q[1];
rz(-1.5412953) q[1];
sx q[1];
rz(-2.704481) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1431883) q[0];
sx q[0];
rz(-2.8371053) q[0];
sx q[0];
rz(1.2782774) q[0];
rz(-pi) q[1];
rz(-2.897981) q[2];
sx q[2];
rz(-1.1738452) q[2];
sx q[2];
rz(1.5543907) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.25698173) q[1];
sx q[1];
rz(-1.4323438) q[1];
sx q[1];
rz(-0.42582663) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0760355) q[3];
sx q[3];
rz(-1.5177739) q[3];
sx q[3];
rz(3.0474636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1221216) q[2];
sx q[2];
rz(-0.78142587) q[2];
sx q[2];
rz(2.1195228) q[2];
rz(1.2381037) q[3];
sx q[3];
rz(-2.759203) q[3];
sx q[3];
rz(-2.7220272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5220752) q[0];
sx q[0];
rz(-1.8958805) q[0];
sx q[0];
rz(2.1602901) q[0];
rz(3.006382) q[1];
sx q[1];
rz(-1.0842666) q[1];
sx q[1];
rz(-2.9503126) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080973074) q[0];
sx q[0];
rz(-1.7195716) q[0];
sx q[0];
rz(3.0540375) q[0];
rz(0.83696604) q[2];
sx q[2];
rz(-1.776473) q[2];
sx q[2];
rz(-2.2107746) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8297255) q[1];
sx q[1];
rz(-0.79677454) q[1];
sx q[1];
rz(-1.7142678) q[1];
rz(-pi) q[2];
rz(-2.5883834) q[3];
sx q[3];
rz(-2.534453) q[3];
sx q[3];
rz(0.34724423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.68025756) q[2];
sx q[2];
rz(-2.1562083) q[2];
sx q[2];
rz(2.130924) q[2];
rz(-2.3800395) q[3];
sx q[3];
rz(-1.9617617) q[3];
sx q[3];
rz(-0.23553577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8028832) q[0];
sx q[0];
rz(-0.25512472) q[0];
sx q[0];
rz(0.55661911) q[0];
rz(3.026475) q[1];
sx q[1];
rz(-1.8042253) q[1];
sx q[1];
rz(-2.1690878) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33362493) q[0];
sx q[0];
rz(-1.4687612) q[0];
sx q[0];
rz(1.5029961) q[0];
rz(-pi) q[1];
rz(-0.82053484) q[2];
sx q[2];
rz(-1.3247326) q[2];
sx q[2];
rz(-2.8503502) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.16469615) q[1];
sx q[1];
rz(-1.4310734) q[1];
sx q[1];
rz(1.7663899) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1049018) q[3];
sx q[3];
rz(-2.1864236) q[3];
sx q[3];
rz(-0.43945593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3187023) q[2];
sx q[2];
rz(-2.0501142) q[2];
sx q[2];
rz(-1.5931607) q[2];
rz(1.7758153) q[3];
sx q[3];
rz(-2.8184991) q[3];
sx q[3];
rz(-2.1877066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4218629) q[0];
sx q[0];
rz(-1.8780163) q[0];
sx q[0];
rz(1.7156037) q[0];
rz(-1.0643719) q[1];
sx q[1];
rz(-2.1247037) q[1];
sx q[1];
rz(2.7672966) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74852809) q[0];
sx q[0];
rz(-0.52951282) q[0];
sx q[0];
rz(-0.66324309) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3099953) q[2];
sx q[2];
rz(-2.0506095) q[2];
sx q[2];
rz(-0.11670437) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4243851) q[1];
sx q[1];
rz(-1.8339515) q[1];
sx q[1];
rz(-0.56916635) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4738594) q[3];
sx q[3];
rz(-2.1453834) q[3];
sx q[3];
rz(2.7001065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8950243) q[2];
sx q[2];
rz(-0.66528577) q[2];
sx q[2];
rz(-2.1833615) q[2];
rz(2.9124177) q[3];
sx q[3];
rz(-1.6848247) q[3];
sx q[3];
rz(-2.5206101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36528698) q[0];
sx q[0];
rz(-1.1927274) q[0];
sx q[0];
rz(0.90674415) q[0];
rz(-1.0892185) q[1];
sx q[1];
rz(-1.6420495) q[1];
sx q[1];
rz(1.8315171) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3416672) q[0];
sx q[0];
rz(-2.7503715) q[0];
sx q[0];
rz(-0.13473405) q[0];
x q[1];
rz(-1.2398948) q[2];
sx q[2];
rz(-0.885193) q[2];
sx q[2];
rz(-2.9943525) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.29855) q[1];
sx q[1];
rz(-1.3805461) q[1];
sx q[1];
rz(-0.33903867) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0974) q[3];
sx q[3];
rz(-2.5857539) q[3];
sx q[3];
rz(1.0122055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.92581302) q[2];
sx q[2];
rz(-0.44162193) q[2];
sx q[2];
rz(1.4833935) q[2];
rz(-0.27967927) q[3];
sx q[3];
rz(-0.99273434) q[3];
sx q[3];
rz(-0.057597615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.9782372) q[0];
sx q[0];
rz(-1.5196479) q[0];
sx q[0];
rz(-0.21959198) q[0];
rz(0.50312463) q[1];
sx q[1];
rz(-0.88880912) q[1];
sx q[1];
rz(-0.84987744) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69707623) q[0];
sx q[0];
rz(-0.27643099) q[0];
sx q[0];
rz(-0.96093775) q[0];
rz(3.0476961) q[2];
sx q[2];
rz(-2.2517423) q[2];
sx q[2];
rz(-0.64881334) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5819468) q[1];
sx q[1];
rz(-1.5362745) q[1];
sx q[1];
rz(0.80969651) q[1];
x q[2];
rz(1.3109342) q[3];
sx q[3];
rz(-1.6901008) q[3];
sx q[3];
rz(0.80727778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5510817) q[2];
sx q[2];
rz(-1.3122281) q[2];
sx q[2];
rz(-1.760651) q[2];
rz(-2.3855709) q[3];
sx q[3];
rz(-0.20320007) q[3];
sx q[3];
rz(0.35593629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6324156) q[0];
sx q[0];
rz(-0.89288765) q[0];
sx q[0];
rz(2.7365622) q[0];
rz(0.45267725) q[1];
sx q[1];
rz(-2.15937) q[1];
sx q[1];
rz(-1.8639494) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4687846) q[0];
sx q[0];
rz(-1.0793669) q[0];
sx q[0];
rz(1.8322893) q[0];
rz(0.69127609) q[2];
sx q[2];
rz(-1.9553767) q[2];
sx q[2];
rz(-0.28011766) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9459878) q[1];
sx q[1];
rz(-2.0161649) q[1];
sx q[1];
rz(-3.1006378) q[1];
rz(-pi) q[2];
rz(-2.3074179) q[3];
sx q[3];
rz(-1.8636384) q[3];
sx q[3];
rz(1.0398231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8032802) q[2];
sx q[2];
rz(-2.7351604) q[2];
sx q[2];
rz(2.7837616) q[2];
rz(-1.7221649) q[3];
sx q[3];
rz(-1.2737041) q[3];
sx q[3];
rz(-2.0675802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91838592) q[0];
sx q[0];
rz(-3.0637488) q[0];
sx q[0];
rz(3.0293368) q[0];
rz(-0.90011251) q[1];
sx q[1];
rz(-1.0670412) q[1];
sx q[1];
rz(0.21044883) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2699379) q[0];
sx q[0];
rz(-2.6007814) q[0];
sx q[0];
rz(0.35738118) q[0];
rz(-pi) q[1];
rz(-1.9824355) q[2];
sx q[2];
rz(-1.4201846) q[2];
sx q[2];
rz(2.6580236) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3810972) q[1];
sx q[1];
rz(-1.3638745) q[1];
sx q[1];
rz(-0.31462545) q[1];
rz(2.3841303) q[3];
sx q[3];
rz(-2.8031581) q[3];
sx q[3];
rz(-2.7528742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9615053) q[2];
sx q[2];
rz(-2.4752361) q[2];
sx q[2];
rz(-1.5853184) q[2];
rz(1.8680343) q[3];
sx q[3];
rz(-0.62265101) q[3];
sx q[3];
rz(-0.20726985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56959854) q[0];
sx q[0];
rz(-0.8710237) q[0];
sx q[0];
rz(-1.3652753) q[0];
rz(2.3251484) q[1];
sx q[1];
rz(-1.2533617) q[1];
sx q[1];
rz(-0.15773699) q[1];
rz(2.9755637) q[2];
sx q[2];
rz(-1.1062853) q[2];
sx q[2];
rz(1.7044978) q[2];
rz(-2.7995085) q[3];
sx q[3];
rz(-0.87089201) q[3];
sx q[3];
rz(1.0637829) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];