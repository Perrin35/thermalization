OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8776956) q[0];
sx q[0];
rz(-2.1030302) q[0];
sx q[0];
rz(-2.7466018) q[0];
rz(0.78139961) q[1];
sx q[1];
rz(-1.7698987) q[1];
sx q[1];
rz(2.3514907) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27777729) q[0];
sx q[0];
rz(-0.57584563) q[0];
sx q[0];
rz(-3.083549) q[0];
rz(-pi) q[1];
rz(1.1320791) q[2];
sx q[2];
rz(-0.84442511) q[2];
sx q[2];
rz(1.7199291) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2826281) q[1];
sx q[1];
rz(-0.55759768) q[1];
sx q[1];
rz(-1.2931812) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4710841) q[3];
sx q[3];
rz(-2.5250348) q[3];
sx q[3];
rz(-2.3634865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6140952) q[2];
sx q[2];
rz(-3.0588394) q[2];
sx q[2];
rz(-2.5772074) q[2];
rz(-2.2218521) q[3];
sx q[3];
rz(-2.4904833) q[3];
sx q[3];
rz(-1.1434327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1564002) q[0];
sx q[0];
rz(-2.1512478) q[0];
sx q[0];
rz(-1.9957969) q[0];
rz(1.0076373) q[1];
sx q[1];
rz(-2.2538908) q[1];
sx q[1];
rz(0.064858286) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.295395) q[0];
sx q[0];
rz(-1.0480651) q[0];
sx q[0];
rz(0.46434648) q[0];
rz(-pi) q[1];
rz(0.95352651) q[2];
sx q[2];
rz(-0.89260403) q[2];
sx q[2];
rz(1.2645406) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8255348) q[1];
sx q[1];
rz(-1.59606) q[1];
sx q[1];
rz(0.18417321) q[1];
rz(-1.0367619) q[3];
sx q[3];
rz(-2.0151842) q[3];
sx q[3];
rz(2.5787974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2603944) q[2];
sx q[2];
rz(-0.10513267) q[2];
sx q[2];
rz(-2.6283188) q[2];
rz(-1.5848292) q[3];
sx q[3];
rz(-1.9624174) q[3];
sx q[3];
rz(-1.56196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8455115) q[0];
sx q[0];
rz(-0.054231461) q[0];
sx q[0];
rz(-0.84501141) q[0];
rz(-2.4301279) q[1];
sx q[1];
rz(-0.80159801) q[1];
sx q[1];
rz(2.5149288) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1820647) q[0];
sx q[0];
rz(-1.7668889) q[0];
sx q[0];
rz(0.0053868731) q[0];
rz(-pi) q[1];
x q[1];
rz(1.955613) q[2];
sx q[2];
rz(-1.3582412) q[2];
sx q[2];
rz(0.40639588) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2528673) q[1];
sx q[1];
rz(-1.3952726) q[1];
sx q[1];
rz(-0.85078001) q[1];
rz(-pi) q[2];
rz(-0.0073284433) q[3];
sx q[3];
rz(-1.0798608) q[3];
sx q[3];
rz(-1.0388663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6381548) q[2];
sx q[2];
rz(-1.902709) q[2];
sx q[2];
rz(2.9729291) q[2];
rz(-2.9851798) q[3];
sx q[3];
rz(-0.63786879) q[3];
sx q[3];
rz(2.8580581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4446568) q[0];
sx q[0];
rz(-2.038027) q[0];
sx q[0];
rz(-2.7705833) q[0];
rz(-1.8095398) q[1];
sx q[1];
rz(-1.5947554) q[1];
sx q[1];
rz(0.36088774) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9103521) q[0];
sx q[0];
rz(-0.12847289) q[0];
sx q[0];
rz(-0.84004004) q[0];
rz(-pi) q[1];
rz(0.12207403) q[2];
sx q[2];
rz(-1.5549107) q[2];
sx q[2];
rz(1.0779013) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0119355) q[1];
sx q[1];
rz(-1.8399554) q[1];
sx q[1];
rz(0.98711826) q[1];
x q[2];
rz(2.1122574) q[3];
sx q[3];
rz(-0.96857054) q[3];
sx q[3];
rz(2.9871939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0767625) q[2];
sx q[2];
rz(-2.7702489) q[2];
sx q[2];
rz(-1.4999464) q[2];
rz(2.0867945) q[3];
sx q[3];
rz(-1.0902371) q[3];
sx q[3];
rz(-1.1928026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.763279) q[0];
sx q[0];
rz(-1.9739456) q[0];
sx q[0];
rz(1.8462697) q[0];
rz(3.0150705) q[1];
sx q[1];
rz(-2.449072) q[1];
sx q[1];
rz(1.893938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0678789) q[0];
sx q[0];
rz(-2.824221) q[0];
sx q[0];
rz(0.48010357) q[0];
x q[1];
rz(-0.7681925) q[2];
sx q[2];
rz(-1.2821226) q[2];
sx q[2];
rz(-1.7199668) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.32198731) q[1];
sx q[1];
rz(-1.7778559) q[1];
sx q[1];
rz(-0.38298573) q[1];
rz(1.2887824) q[3];
sx q[3];
rz(-1.7502893) q[3];
sx q[3];
rz(3.0075551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7724472) q[2];
sx q[2];
rz(-2.7771066) q[2];
sx q[2];
rz(-0.44949731) q[2];
rz(0.43826023) q[3];
sx q[3];
rz(-2.0354383) q[3];
sx q[3];
rz(2.0800169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19675572) q[0];
sx q[0];
rz(-1.0642835) q[0];
sx q[0];
rz(0.34360316) q[0];
rz(-2.4876439) q[1];
sx q[1];
rz(-2.3908354) q[1];
sx q[1];
rz(-2.5508945) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66309975) q[0];
sx q[0];
rz(-1.0859153) q[0];
sx q[0];
rz(-0.25867489) q[0];
x q[1];
rz(-0.74004897) q[2];
sx q[2];
rz(-2.165926) q[2];
sx q[2];
rz(2.8071432) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0456344) q[1];
sx q[1];
rz(-1.5603298) q[1];
sx q[1];
rz(-1.7016618) q[1];
rz(-2.2600987) q[3];
sx q[3];
rz(-2.1432264) q[3];
sx q[3];
rz(1.5015757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.58064738) q[2];
sx q[2];
rz(-0.86981589) q[2];
sx q[2];
rz(-0.84442863) q[2];
rz(-2.1010418) q[3];
sx q[3];
rz(-1.8620164) q[3];
sx q[3];
rz(-1.9677264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3867253) q[0];
sx q[0];
rz(-2.6193021) q[0];
sx q[0];
rz(-2.0544384) q[0];
rz(-2.6548751) q[1];
sx q[1];
rz(-1.4954647) q[1];
sx q[1];
rz(-1.6513991) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54571153) q[0];
sx q[0];
rz(-2.427248) q[0];
sx q[0];
rz(1.7313135) q[0];
x q[1];
rz(1.0045687) q[2];
sx q[2];
rz(-0.90551584) q[2];
sx q[2];
rz(-0.69349223) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.76032066) q[1];
sx q[1];
rz(-0.57653344) q[1];
sx q[1];
rz(-0.74200304) q[1];
rz(-2.1726726) q[3];
sx q[3];
rz(-0.35588297) q[3];
sx q[3];
rz(2.8371135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3460441) q[2];
sx q[2];
rz(-2.2717924) q[2];
sx q[2];
rz(2.8022433) q[2];
rz(0.36335534) q[3];
sx q[3];
rz(-0.27000579) q[3];
sx q[3];
rz(1.5800765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5245847) q[0];
sx q[0];
rz(-2.0211077) q[0];
sx q[0];
rz(0.28550112) q[0];
rz(-0.78954804) q[1];
sx q[1];
rz(-1.5317081) q[1];
sx q[1];
rz(-1.5910566) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94211468) q[0];
sx q[0];
rz(-1.7492332) q[0];
sx q[0];
rz(0.30233827) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8218846) q[2];
sx q[2];
rz(-2.6007915) q[2];
sx q[2];
rz(-0.021440949) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4414822) q[1];
sx q[1];
rz(-0.69789819) q[1];
sx q[1];
rz(-0.68241193) q[1];
x q[2];
rz(3.0394394) q[3];
sx q[3];
rz(-2.4718067) q[3];
sx q[3];
rz(2.6492491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.32621128) q[2];
sx q[2];
rz(-0.15924328) q[2];
sx q[2];
rz(2.2747269) q[2];
rz(-2.4032812) q[3];
sx q[3];
rz(-1.5840931) q[3];
sx q[3];
rz(-2.2332634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065780491) q[0];
sx q[0];
rz(-0.74956885) q[0];
sx q[0];
rz(2.442389) q[0];
rz(0.54222822) q[1];
sx q[1];
rz(-1.3424073) q[1];
sx q[1];
rz(0.29621616) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86481536) q[0];
sx q[0];
rz(-0.84057284) q[0];
sx q[0];
rz(0.093412799) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3260957) q[2];
sx q[2];
rz(-1.0692714) q[2];
sx q[2];
rz(-3.0072711) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6275639) q[1];
sx q[1];
rz(-0.96935287) q[1];
sx q[1];
rz(-0.87694962) q[1];
x q[2];
rz(-0.086478905) q[3];
sx q[3];
rz(-1.9473377) q[3];
sx q[3];
rz(1.2264156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2518623) q[2];
sx q[2];
rz(-0.44318649) q[2];
sx q[2];
rz(2.6207793) q[2];
rz(3.1098747) q[3];
sx q[3];
rz(-0.40516502) q[3];
sx q[3];
rz(3.090455) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6601335) q[0];
sx q[0];
rz(-1.0986468) q[0];
sx q[0];
rz(0.81137401) q[0];
rz(-0.42586455) q[1];
sx q[1];
rz(-0.90161294) q[1];
sx q[1];
rz(-1.3070377) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9872914) q[0];
sx q[0];
rz(-1.5736921) q[0];
sx q[0];
rz(-1.5737246) q[0];
rz(2.0352896) q[2];
sx q[2];
rz(-1.5519605) q[2];
sx q[2];
rz(-1.0412168) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9237808) q[1];
sx q[1];
rz(-0.56107035) q[1];
sx q[1];
rz(-0.63528676) q[1];
x q[2];
rz(1.3974992) q[3];
sx q[3];
rz(-1.2504645) q[3];
sx q[3];
rz(0.082435247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1618774) q[2];
sx q[2];
rz(-0.76354176) q[2];
sx q[2];
rz(0.36583501) q[2];
rz(-2.2063935) q[3];
sx q[3];
rz(-1.5465982) q[3];
sx q[3];
rz(2.759815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81015051) q[0];
sx q[0];
rz(-1.3642385) q[0];
sx q[0];
rz(1.4862899) q[0];
rz(1.4797795) q[1];
sx q[1];
rz(-1.231709) q[1];
sx q[1];
rz(-0.92373893) q[1];
rz(0.059546373) q[2];
sx q[2];
rz(-2.1179143) q[2];
sx q[2];
rz(2.2512022) q[2];
rz(-2.6213358) q[3];
sx q[3];
rz(-2.8388966) q[3];
sx q[3];
rz(2.7404529) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
