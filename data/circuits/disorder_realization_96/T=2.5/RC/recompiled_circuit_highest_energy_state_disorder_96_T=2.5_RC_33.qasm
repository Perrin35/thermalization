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
rz(1.9198298) q[0];
sx q[0];
rz(-1.3718995) q[0];
sx q[0];
rz(-2.0576117) q[0];
rz(2.3594175) q[1];
sx q[1];
rz(-0.44936925) q[1];
sx q[1];
rz(2.3545797) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5840658) q[0];
sx q[0];
rz(-1.9587933) q[0];
sx q[0];
rz(-0.16325133) q[0];
x q[1];
rz(2.2465503) q[2];
sx q[2];
rz(-1.6056201) q[2];
sx q[2];
rz(0.36961242) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1842332) q[1];
sx q[1];
rz(-2.5788947) q[1];
sx q[1];
rz(1.301728) q[1];
rz(-pi) q[2];
x q[2];
rz(1.161133) q[3];
sx q[3];
rz(-0.72736638) q[3];
sx q[3];
rz(-1.1687708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3419753) q[2];
sx q[2];
rz(-1.3222597) q[2];
sx q[2];
rz(-0.70023099) q[2];
rz(1.7796984) q[3];
sx q[3];
rz(-1.9340065) q[3];
sx q[3];
rz(-0.12968682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.48548651) q[0];
sx q[0];
rz(-2.7950197) q[0];
sx q[0];
rz(-1.19278) q[0];
rz(2.9876409) q[1];
sx q[1];
rz(-0.62869453) q[1];
sx q[1];
rz(-1.2492294) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54493413) q[0];
sx q[0];
rz(-1.3976263) q[0];
sx q[0];
rz(1.5073266) q[0];
x q[1];
rz(-3.0871816) q[2];
sx q[2];
rz(-0.93428946) q[2];
sx q[2];
rz(-2.7639703) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.8802048) q[1];
sx q[1];
rz(-1.4900165) q[1];
sx q[1];
rz(-0.30552633) q[1];
rz(0.90460299) q[3];
sx q[3];
rz(-1.8663916) q[3];
sx q[3];
rz(-1.4669661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.93231702) q[2];
sx q[2];
rz(-1.2195769) q[2];
sx q[2];
rz(-0.90927124) q[2];
rz(-3.0026109) q[3];
sx q[3];
rz(-0.2551955) q[3];
sx q[3];
rz(-2.2709258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9210377) q[0];
sx q[0];
rz(-1.9746566) q[0];
sx q[0];
rz(-1.1880818) q[0];
rz(2.2782169) q[1];
sx q[1];
rz(-1.2438351) q[1];
sx q[1];
rz(-0.9695425) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.506451) q[0];
sx q[0];
rz(-1.4124845) q[0];
sx q[0];
rz(-0.16198762) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9447504) q[2];
sx q[2];
rz(-1.2450758) q[2];
sx q[2];
rz(-1.0870799) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0672966) q[1];
sx q[1];
rz(-2.887425) q[1];
sx q[1];
rz(-0.39293082) q[1];
rz(-2.5695877) q[3];
sx q[3];
rz(-2.1915132) q[3];
sx q[3];
rz(-2.2748924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.36677507) q[2];
sx q[2];
rz(-0.25664169) q[2];
sx q[2];
rz(-0.65566629) q[2];
rz(1.5244779) q[3];
sx q[3];
rz(-1.7888125) q[3];
sx q[3];
rz(-0.85339439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26576385) q[0];
sx q[0];
rz(-2.0142856) q[0];
sx q[0];
rz(1.458459) q[0];
rz(1.726285) q[1];
sx q[1];
rz(-1.7067319) q[1];
sx q[1];
rz(1.2687792) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44132933) q[0];
sx q[0];
rz(-0.95281314) q[0];
sx q[0];
rz(2.000314) q[0];
x q[1];
rz(0.94082035) q[2];
sx q[2];
rz(-2.0220304) q[2];
sx q[2];
rz(-2.4521417) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3038588) q[1];
sx q[1];
rz(-1.7540252) q[1];
sx q[1];
rz(0.86838037) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7990718) q[3];
sx q[3];
rz(-1.6581495) q[3];
sx q[3];
rz(-0.5783847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.65712523) q[2];
sx q[2];
rz(-1.7871658) q[2];
sx q[2];
rz(-1.8787059) q[2];
rz(-0.2028939) q[3];
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
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6577067) q[0];
sx q[0];
rz(-1.0840451) q[0];
sx q[0];
rz(-0.44802353) q[0];
rz(-1.1605284) q[1];
sx q[1];
rz(-1.0803761) q[1];
sx q[1];
rz(-2.0652658) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1370476) q[0];
sx q[0];
rz(-2.2820615) q[0];
sx q[0];
rz(0.28202951) q[0];
rz(0.14400173) q[2];
sx q[2];
rz(-0.33641923) q[2];
sx q[2];
rz(-0.1937723) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4566735) q[1];
sx q[1];
rz(-3.0701048) q[1];
sx q[1];
rz(2.100657) q[1];
rz(-pi) q[2];
rz(-0.95376261) q[3];
sx q[3];
rz(-1.5542398) q[3];
sx q[3];
rz(-1.7615033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6059604) q[2];
sx q[2];
rz(-1.2781906) q[2];
sx q[2];
rz(-0.34206259) q[2];
rz(-2.0096807) q[3];
sx q[3];
rz(-1.071238) q[3];
sx q[3];
rz(2.7560077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015942052) q[0];
sx q[0];
rz(-0.60369879) q[0];
sx q[0];
rz(2.0397256) q[0];
rz(-2.8562538) q[1];
sx q[1];
rz(-2.1632748) q[1];
sx q[1];
rz(-2.2314821) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21134899) q[0];
sx q[0];
rz(-1.4574629) q[0];
sx q[0];
rz(1.428753) q[0];
x q[1];
rz(3.0481045) q[2];
sx q[2];
rz(-1.4564691) q[2];
sx q[2];
rz(-1.0766407) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6807947) q[1];
sx q[1];
rz(-2.6943047) q[1];
sx q[1];
rz(-1.4385002) q[1];
x q[2];
rz(1.7496787) q[3];
sx q[3];
rz(-1.4174882) q[3];
sx q[3];
rz(-1.2962504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.82112271) q[2];
sx q[2];
rz(-1.4218825) q[2];
sx q[2];
rz(-1.0538496) q[2];
rz(-0.36367917) q[3];
sx q[3];
rz(-1.6561008) q[3];
sx q[3];
rz(0.53500879) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66154552) q[0];
sx q[0];
rz(-2.1286025) q[0];
sx q[0];
rz(-2.9841828) q[0];
rz(-2.5120381) q[1];
sx q[1];
rz(-1.5512356) q[1];
sx q[1];
rz(-3.0622283) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.019477) q[0];
sx q[0];
rz(-1.6412897) q[0];
sx q[0];
rz(1.6623783) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3936226) q[2];
sx q[2];
rz(-1.9594155) q[2];
sx q[2];
rz(1.7541698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.094202894) q[1];
sx q[1];
rz(-1.4416639) q[1];
sx q[1];
rz(-2.6455621) q[1];
x q[2];
rz(0.98439545) q[3];
sx q[3];
rz(-1.5933872) q[3];
sx q[3];
rz(0.65832321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6307512) q[2];
sx q[2];
rz(-1.8264822) q[2];
sx q[2];
rz(-2.9691479) q[2];
rz(2.4991961) q[3];
sx q[3];
rz(-2.1746217) q[3];
sx q[3];
rz(0.35014686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
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
rz(-2.9134275) q[1];
sx q[1];
rz(-1.4527495) q[1];
sx q[1];
rz(-1.3120922) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23240664) q[0];
sx q[0];
rz(-2.3433422) q[0];
sx q[0];
rz(-1.2212103) q[0];
x q[1];
rz(1.4390575) q[2];
sx q[2];
rz(-2.4102168) q[2];
sx q[2];
rz(0.057550427) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6375288) q[1];
sx q[1];
rz(-0.78128587) q[1];
sx q[1];
rz(0.81721304) q[1];
x q[2];
rz(-0.033574738) q[3];
sx q[3];
rz(-1.9596425) q[3];
sx q[3];
rz(-1.7186708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.18099004) q[2];
sx q[2];
rz(-1.484551) q[2];
sx q[2];
rz(0.20299882) q[2];
rz(2.3024043) q[3];
sx q[3];
rz(-2.8036717) q[3];
sx q[3];
rz(-0.85116974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39076552) q[0];
sx q[0];
rz(-0.44742328) q[0];
sx q[0];
rz(-0.14697337) q[0];
rz(2.7087063) q[1];
sx q[1];
rz(-1.1914445) q[1];
sx q[1];
rz(-1.8870707) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0085274335) q[0];
sx q[0];
rz(-2.0957895) q[0];
sx q[0];
rz(-1.9740482) q[0];
x q[1];
rz(-0.033795088) q[2];
sx q[2];
rz(-1.8189478) q[2];
sx q[2];
rz(-0.42873086) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7045198) q[1];
sx q[1];
rz(-2.7531248) q[1];
sx q[1];
rz(1.1004431) q[1];
rz(-1.5490554) q[3];
sx q[3];
rz(-1.4026902) q[3];
sx q[3];
rz(-0.67847195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7811084) q[2];
sx q[2];
rz(-1.0177178) q[2];
sx q[2];
rz(2.600889) q[2];
rz(-0.30132076) q[3];
sx q[3];
rz(-0.40073985) q[3];
sx q[3];
rz(-1.8436684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8526469) q[0];
sx q[0];
rz(-1.1971373) q[0];
sx q[0];
rz(2.4135015) q[0];
rz(2.6241265) q[1];
sx q[1];
rz(-1.6518075) q[1];
sx q[1];
rz(-0.99658406) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9835127) q[0];
sx q[0];
rz(-1.5444813) q[0];
sx q[0];
rz(1.2264432) q[0];
rz(-pi) q[1];
rz(0.58169843) q[2];
sx q[2];
rz(-1.4862183) q[2];
sx q[2];
rz(-2.4385114) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7000632) q[1];
sx q[1];
rz(-0.25694381) q[1];
sx q[1];
rz(-1.5154626) q[1];
rz(-pi) q[2];
rz(-0.49788614) q[3];
sx q[3];
rz(-0.20789805) q[3];
sx q[3];
rz(2.9836536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1835798) q[2];
sx q[2];
rz(-1.5637584) q[2];
sx q[2];
rz(-3.0713522) q[2];
rz(3.1349414) q[3];
sx q[3];
rz(-2.9700322) q[3];
sx q[3];
rz(-3.0370965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9329279) q[0];
sx q[0];
rz(-1.4069955) q[0];
sx q[0];
rz(1.7250742) q[0];
rz(0.76864645) q[1];
sx q[1];
rz(-1.2223163) q[1];
sx q[1];
rz(-0.32739122) q[1];
rz(-3.1216172) q[2];
sx q[2];
rz(-1.6835313) q[2];
sx q[2];
rz(2.6575553) q[2];
rz(-1.205659) q[3];
sx q[3];
rz(-1.3177928) q[3];
sx q[3];
rz(-3.1402088) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
