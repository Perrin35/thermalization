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
rz(-1.2217628) q[0];
sx q[0];
rz(-1.7696932) q[0];
sx q[0];
rz(2.0576117) q[0];
rz(2.3594175) q[1];
sx q[1];
rz(-0.44936925) q[1];
sx q[1];
rz(2.3545797) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9946905) q[0];
sx q[0];
rz(-2.7222607) q[0];
sx q[0];
rz(-1.9492784) q[0];
rz(1.6264362) q[2];
sx q[2];
rz(-0.67650992) q[2];
sx q[2];
rz(-1.8970053) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38429425) q[1];
sx q[1];
rz(-1.7130902) q[1];
sx q[1];
rz(1.0245002) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9804597) q[3];
sx q[3];
rz(-0.72736638) q[3];
sx q[3];
rz(-1.9728219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7996173) q[2];
sx q[2];
rz(-1.3222597) q[2];
sx q[2];
rz(2.4413617) q[2];
rz(-1.7796984) q[3];
sx q[3];
rz(-1.9340065) q[3];
sx q[3];
rz(-3.0119058) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48548651) q[0];
sx q[0];
rz(-0.34657297) q[0];
sx q[0];
rz(1.19278) q[0];
rz(-2.9876409) q[1];
sx q[1];
rz(-2.5128981) q[1];
sx q[1];
rz(1.8923632) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0368127) q[0];
sx q[0];
rz(-1.6333155) q[0];
sx q[0];
rz(0.17351242) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0871816) q[2];
sx q[2];
rz(-2.2073032) q[2];
sx q[2];
rz(-2.7639703) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4255562) q[1];
sx q[1];
rz(-1.2662983) q[1];
sx q[1];
rz(-1.4861121) q[1];
x q[2];
rz(2.2369897) q[3];
sx q[3];
rz(-1.8663916) q[3];
sx q[3];
rz(1.4669661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93231702) q[2];
sx q[2];
rz(-1.2195769) q[2];
sx q[2];
rz(-0.90927124) q[2];
rz(-0.13898177) q[3];
sx q[3];
rz(-2.8863972) q[3];
sx q[3];
rz(0.87066686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22055498) q[0];
sx q[0];
rz(-1.9746566) q[0];
sx q[0];
rz(1.9535109) q[0];
rz(-0.86337572) q[1];
sx q[1];
rz(-1.2438351) q[1];
sx q[1];
rz(2.1720502) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6351417) q[0];
sx q[0];
rz(-1.4124845) q[0];
sx q[0];
rz(-0.16198762) q[0];
rz(-pi) q[1];
rz(-1.1968422) q[2];
sx q[2];
rz(-1.2450758) q[2];
sx q[2];
rz(-2.0545127) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.66268209) q[1];
sx q[1];
rz(-1.3363774) q[1];
sx q[1];
rz(-1.6699416) q[1];
rz(-pi) q[2];
rz(0.86609716) q[3];
sx q[3];
rz(-2.0267762) q[3];
sx q[3];
rz(2.7957819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.36677507) q[2];
sx q[2];
rz(-0.25664169) q[2];
sx q[2];
rz(0.65566629) q[2];
rz(1.5244779) q[3];
sx q[3];
rz(-1.7888125) q[3];
sx q[3];
rz(-0.85339439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26576385) q[0];
sx q[0];
rz(-2.0142856) q[0];
sx q[0];
rz(1.458459) q[0];
rz(1.4153076) q[1];
sx q[1];
rz(-1.4348607) q[1];
sx q[1];
rz(1.2687792) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7527217) q[0];
sx q[0];
rz(-1.9170887) q[0];
sx q[0];
rz(-2.4779921) q[0];
rz(0.54019467) q[2];
sx q[2];
rz(-2.1295851) q[2];
sx q[2];
rz(1.9524198) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8377338) q[1];
sx q[1];
rz(-1.3875675) q[1];
sx q[1];
rz(-2.2732123) q[1];
rz(-pi) q[2];
rz(-1.7990718) q[3];
sx q[3];
rz(-1.4834431) q[3];
sx q[3];
rz(2.563208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4844674) q[2];
sx q[2];
rz(-1.3544269) q[2];
sx q[2];
rz(-1.8787059) q[2];
rz(2.9386988) q[3];
sx q[3];
rz(-1.032136) q[3];
sx q[3];
rz(-1.3304905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.48388594) q[0];
sx q[0];
rz(-2.0575476) q[0];
sx q[0];
rz(-2.6935691) q[0];
rz(-1.1605284) q[1];
sx q[1];
rz(-2.0612165) q[1];
sx q[1];
rz(-1.0763268) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1370476) q[0];
sx q[0];
rz(-0.8595312) q[0];
sx q[0];
rz(2.8595631) q[0];
rz(-pi) q[1];
rz(0.14400173) q[2];
sx q[2];
rz(-2.8051734) q[2];
sx q[2];
rz(0.1937723) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68491918) q[1];
sx q[1];
rz(-3.0701048) q[1];
sx q[1];
rz(-1.0409357) q[1];
rz(-pi) q[2];
rz(-0.020298941) q[3];
sx q[3];
rz(-0.95385984) q[3];
sx q[3];
rz(2.9391409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5356323) q[2];
sx q[2];
rz(-1.863402) q[2];
sx q[2];
rz(-0.34206259) q[2];
rz(-2.0096807) q[3];
sx q[3];
rz(-2.0703546) q[3];
sx q[3];
rz(-2.7560077) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015942052) q[0];
sx q[0];
rz(-0.60369879) q[0];
sx q[0];
rz(-1.1018671) q[0];
rz(0.28533882) q[1];
sx q[1];
rz(-2.1632748) q[1];
sx q[1];
rz(0.91011059) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21134899) q[0];
sx q[0];
rz(-1.6841297) q[0];
sx q[0];
rz(1.7128397) q[0];
rz(-pi) q[1];
rz(-0.88819529) q[2];
sx q[2];
rz(-2.9940372) q[2];
sx q[2];
rz(0.38868586) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.460798) q[1];
sx q[1];
rz(-2.6943047) q[1];
sx q[1];
rz(-1.4385002) q[1];
rz(-0.85569622) q[3];
sx q[3];
rz(-0.2350546) q[3];
sx q[3];
rz(-2.1657221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3204699) q[2];
sx q[2];
rz(-1.7197101) q[2];
sx q[2];
rz(-2.087743) q[2];
rz(2.7779135) q[3];
sx q[3];
rz(-1.6561008) q[3];
sx q[3];
rz(-2.6065839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66154552) q[0];
sx q[0];
rz(-1.0129901) q[0];
sx q[0];
rz(-0.15740982) q[0];
rz(2.5120381) q[1];
sx q[1];
rz(-1.5512356) q[1];
sx q[1];
rz(-0.079364337) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12211563) q[0];
sx q[0];
rz(-1.6412897) q[0];
sx q[0];
rz(-1.4792144) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3936226) q[2];
sx q[2];
rz(-1.9594155) q[2];
sx q[2];
rz(1.7541698) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5461681) q[1];
sx q[1];
rz(-1.079266) q[1];
sx q[1];
rz(-1.4242055) q[1];
rz(-2.1571972) q[3];
sx q[3];
rz(-1.5482054) q[3];
sx q[3];
rz(2.4832694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.51084149) q[2];
sx q[2];
rz(-1.8264822) q[2];
sx q[2];
rz(-2.9691479) q[2];
rz(0.6423966) q[3];
sx q[3];
rz(-2.1746217) q[3];
sx q[3];
rz(-0.35014686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.073535) q[0];
sx q[0];
rz(-1.8386766) q[0];
sx q[0];
rz(-2.5750343) q[0];
rz(-0.22816518) q[1];
sx q[1];
rz(-1.4527495) q[1];
sx q[1];
rz(1.3120922) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5540422) q[0];
sx q[0];
rz(-1.8186079) q[0];
sx q[0];
rz(-2.3378563) q[0];
x q[1];
rz(-3.0242537) q[2];
sx q[2];
rz(-0.84716958) q[2];
sx q[2];
rz(-3.0229508) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6375288) q[1];
sx q[1];
rz(-2.3603068) q[1];
sx q[1];
rz(-0.81721304) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4890461) q[3];
sx q[3];
rz(-0.39021947) q[3];
sx q[3];
rz(-1.8070328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9606026) q[2];
sx q[2];
rz(-1.484551) q[2];
sx q[2];
rz(-0.20299882) q[2];
rz(-0.8391884) q[3];
sx q[3];
rz(-0.33792096) q[3];
sx q[3];
rz(0.85116974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39076552) q[0];
sx q[0];
rz(-2.6941694) q[0];
sx q[0];
rz(-2.9946193) q[0];
rz(-2.7087063) q[1];
sx q[1];
rz(-1.9501481) q[1];
sx q[1];
rz(-1.8870707) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1330652) q[0];
sx q[0];
rz(-1.0458032) q[0];
sx q[0];
rz(-1.9740482) q[0];
x q[1];
rz(-0.033795088) q[2];
sx q[2];
rz(-1.8189478) q[2];
sx q[2];
rz(-0.42873086) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93941037) q[1];
sx q[1];
rz(-1.9152043) q[1];
sx q[1];
rz(-2.9581974) q[1];
rz(-pi) q[2];
rz(-2.9734475) q[3];
sx q[3];
rz(-1.5922308) q[3];
sx q[3];
rz(0.88868621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7811084) q[2];
sx q[2];
rz(-2.1238748) q[2];
sx q[2];
rz(0.54070365) q[2];
rz(-0.30132076) q[3];
sx q[3];
rz(-2.7408528) q[3];
sx q[3];
rz(1.8436684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.28894579) q[0];
sx q[0];
rz(-1.1971373) q[0];
sx q[0];
rz(-0.72809118) q[0];
rz(2.6241265) q[1];
sx q[1];
rz(-1.4897852) q[1];
sx q[1];
rz(0.99658406) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1580799) q[0];
sx q[0];
rz(-1.5971113) q[0];
sx q[0];
rz(1.2264432) q[0];
rz(-1.4696737) q[2];
sx q[2];
rz(-0.99144672) q[2];
sx q[2];
rz(0.92321009) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4987373) q[1];
sx q[1];
rz(-1.3142546) q[1];
sx q[1];
rz(-0.014530583) q[1];
rz(-pi) q[2];
rz(2.9583361) q[3];
sx q[3];
rz(-1.4720639) q[3];
sx q[3];
rz(1.9016641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.95801282) q[2];
sx q[2];
rz(-1.5637584) q[2];
sx q[2];
rz(-0.070240423) q[2];
rz(3.1349414) q[3];
sx q[3];
rz(-0.17156048) q[3];
sx q[3];
rz(3.0370965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2086647) q[0];
sx q[0];
rz(-1.4069955) q[0];
sx q[0];
rz(1.7250742) q[0];
rz(2.3729462) q[1];
sx q[1];
rz(-1.9192764) q[1];
sx q[1];
rz(2.8142014) q[1];
rz(-3.1216172) q[2];
sx q[2];
rz(-1.6835313) q[2];
sx q[2];
rz(2.6575553) q[2];
rz(-0.27003084) q[3];
sx q[3];
rz(-1.217801) q[3];
sx q[3];
rz(1.4767811) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
