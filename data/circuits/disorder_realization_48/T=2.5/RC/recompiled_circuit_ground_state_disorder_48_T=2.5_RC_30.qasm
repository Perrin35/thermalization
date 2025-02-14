OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0410864) q[0];
sx q[0];
rz(-0.20809986) q[0];
sx q[0];
rz(2.7383374) q[0];
rz(2.9085605) q[1];
sx q[1];
rz(-1.7014528) q[1];
sx q[1];
rz(-0.22388248) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96432811) q[0];
sx q[0];
rz(-1.3601662) q[0];
sx q[0];
rz(-0.62646477) q[0];
rz(-pi) q[1];
rz(-3.0133648) q[2];
sx q[2];
rz(-2.2710481) q[2];
sx q[2];
rz(0.70125225) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.76624291) q[1];
sx q[1];
rz(-0.92449429) q[1];
sx q[1];
rz(0.28278858) q[1];
rz(-pi) q[2];
rz(2.3305064) q[3];
sx q[3];
rz(-1.038365) q[3];
sx q[3];
rz(-2.6528751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.69748059) q[2];
sx q[2];
rz(-0.31284249) q[2];
sx q[2];
rz(3.1057788) q[2];
rz(2.4114285) q[3];
sx q[3];
rz(-1.9883479) q[3];
sx q[3];
rz(2.9353976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8908454) q[0];
sx q[0];
rz(-0.99531168) q[0];
sx q[0];
rz(-2.9440951) q[0];
rz(-2.6644871) q[1];
sx q[1];
rz(-0.66480607) q[1];
sx q[1];
rz(0.8055996) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8957414) q[0];
sx q[0];
rz(-0.68994683) q[0];
sx q[0];
rz(0.2941546) q[0];
x q[1];
rz(0.38061541) q[2];
sx q[2];
rz(-1.632649) q[2];
sx q[2];
rz(0.46910367) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14368421) q[1];
sx q[1];
rz(-0.98947462) q[1];
sx q[1];
rz(-2.3475926) q[1];
rz(-pi) q[2];
rz(-2.9904705) q[3];
sx q[3];
rz(-1.7100705) q[3];
sx q[3];
rz(-2.3230011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8311367) q[2];
sx q[2];
rz(-1.6418991) q[2];
sx q[2];
rz(-1.4031225) q[2];
rz(-2.5095615) q[3];
sx q[3];
rz(-2.8681614) q[3];
sx q[3];
rz(-3.0145751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8359351) q[0];
sx q[0];
rz(-2.5129565) q[0];
sx q[0];
rz(-2.054731) q[0];
rz(-0.19733363) q[1];
sx q[1];
rz(-1.5060164) q[1];
sx q[1];
rz(2.8089583) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5088288) q[0];
sx q[0];
rz(-0.45464215) q[0];
sx q[0];
rz(1.2484545) q[0];
x q[1];
rz(-0.83586043) q[2];
sx q[2];
rz(-1.5243013) q[2];
sx q[2];
rz(-2.495386) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.91326317) q[1];
sx q[1];
rz(-0.28013924) q[1];
sx q[1];
rz(2.254451) q[1];
rz(1.8007711) q[3];
sx q[3];
rz(-2.284433) q[3];
sx q[3];
rz(2.7909887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7670333) q[2];
sx q[2];
rz(-1.5908396) q[2];
sx q[2];
rz(1.6732875) q[2];
rz(3.1324006) q[3];
sx q[3];
rz(-0.46949783) q[3];
sx q[3];
rz(-2.3495242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(0.62995768) q[0];
sx q[0];
rz(-0.63794962) q[0];
sx q[0];
rz(1.6988423) q[0];
rz(1.0244145) q[1];
sx q[1];
rz(-1.9099648) q[1];
sx q[1];
rz(-2.3146497) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31931811) q[0];
sx q[0];
rz(-1.7546904) q[0];
sx q[0];
rz(0.075550373) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88281544) q[2];
sx q[2];
rz(-1.4512296) q[2];
sx q[2];
rz(-1.9012698) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0322898) q[1];
sx q[1];
rz(-2.656611) q[1];
sx q[1];
rz(0.46320199) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0204861) q[3];
sx q[3];
rz(-2.2495329) q[3];
sx q[3];
rz(2.8687182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8168762) q[2];
sx q[2];
rz(-0.78410316) q[2];
sx q[2];
rz(1.8640222) q[2];
rz(0.68814021) q[3];
sx q[3];
rz(-2.1161931) q[3];
sx q[3];
rz(-0.96243206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0012896) q[0];
sx q[0];
rz(-1.7004509) q[0];
sx q[0];
rz(-2.9241614) q[0];
rz(-2.8561719) q[1];
sx q[1];
rz(-2.5358584) q[1];
sx q[1];
rz(2.6434456) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4460627) q[0];
sx q[0];
rz(-2.208606) q[0];
sx q[0];
rz(-2.7317156) q[0];
rz(-pi) q[1];
rz(-0.26136036) q[2];
sx q[2];
rz(-1.0685421) q[2];
sx q[2];
rz(-2.4136191) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0951257) q[1];
sx q[1];
rz(-2.6066463) q[1];
sx q[1];
rz(1.1535658) q[1];
rz(-2.7015649) q[3];
sx q[3];
rz(-1.8002568) q[3];
sx q[3];
rz(0.10073951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.51776) q[2];
sx q[2];
rz(-2.0443003) q[2];
sx q[2];
rz(2.9916054) q[2];
rz(-0.013966694) q[3];
sx q[3];
rz(-1.59168) q[3];
sx q[3];
rz(-2.6278031) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6496395) q[0];
sx q[0];
rz(-2.5998901) q[0];
sx q[0];
rz(2.3288222) q[0];
rz(-1.4179519) q[1];
sx q[1];
rz(-1.6287454) q[1];
sx q[1];
rz(-2.0779804) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15062885) q[0];
sx q[0];
rz(-2.0888939) q[0];
sx q[0];
rz(-1.9512779) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16367775) q[2];
sx q[2];
rz(-1.5541758) q[2];
sx q[2];
rz(-1.7134242) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.34142445) q[1];
sx q[1];
rz(-1.582195) q[1];
sx q[1];
rz(-0.59873786) q[1];
rz(2.205414) q[3];
sx q[3];
rz(-0.31226102) q[3];
sx q[3];
rz(-2.0978417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.873988) q[2];
sx q[2];
rz(-2.5817817) q[2];
sx q[2];
rz(1.7248636) q[2];
rz(0.060976107) q[3];
sx q[3];
rz(-0.91106001) q[3];
sx q[3];
rz(1.6772259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7199719) q[0];
sx q[0];
rz(-0.82141972) q[0];
sx q[0];
rz(-0.11909568) q[0];
rz(-0.47856092) q[1];
sx q[1];
rz(-2.3275972) q[1];
sx q[1];
rz(2.4518769) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0942587) q[0];
sx q[0];
rz(-2.2660997) q[0];
sx q[0];
rz(-2.354485) q[0];
x q[1];
rz(-0.48819112) q[2];
sx q[2];
rz(-1.2400157) q[2];
sx q[2];
rz(-2.4622126) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3773681) q[1];
sx q[1];
rz(-0.86478327) q[1];
sx q[1];
rz(3.0824667) q[1];
x q[2];
rz(2.9940412) q[3];
sx q[3];
rz(-2.2780212) q[3];
sx q[3];
rz(0.13015166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0927642) q[2];
sx q[2];
rz(-0.060828716) q[2];
sx q[2];
rz(-1.0678585) q[2];
rz(2.974406) q[3];
sx q[3];
rz(-1.7696295) q[3];
sx q[3];
rz(3.0021477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6086513) q[0];
sx q[0];
rz(-0.81357384) q[0];
sx q[0];
rz(-2.2331878) q[0];
rz(0.99336973) q[1];
sx q[1];
rz(-0.16255957) q[1];
sx q[1];
rz(-0.87432528) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4607166) q[0];
sx q[0];
rz(-1.5881638) q[0];
sx q[0];
rz(3.0094023) q[0];
x q[1];
rz(-1.8285334) q[2];
sx q[2];
rz(-1.589587) q[2];
sx q[2];
rz(-2.7692139) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.711124) q[1];
sx q[1];
rz(-1.5114771) q[1];
sx q[1];
rz(-1.7064246) q[1];
rz(-pi) q[2];
rz(-2.1217974) q[3];
sx q[3];
rz(-1.8809109) q[3];
sx q[3];
rz(1.9581025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.61451644) q[2];
sx q[2];
rz(-0.15915844) q[2];
sx q[2];
rz(0.85421526) q[2];
rz(-2.6863875) q[3];
sx q[3];
rz(-0.56536094) q[3];
sx q[3];
rz(0.2555041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47827569) q[0];
sx q[0];
rz(-1.217696) q[0];
sx q[0];
rz(1.664337) q[0];
rz(1.5669589) q[1];
sx q[1];
rz(-2.4576371) q[1];
sx q[1];
rz(2.4050567) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7304717) q[0];
sx q[0];
rz(-0.53148848) q[0];
sx q[0];
rz(2.0474252) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6746503) q[2];
sx q[2];
rz(-1.1919293) q[2];
sx q[2];
rz(1.9571638) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.354035) q[1];
sx q[1];
rz(-1.327146) q[1];
sx q[1];
rz(0.32275782) q[1];
x q[2];
rz(1.6227342) q[3];
sx q[3];
rz(-0.59762663) q[3];
sx q[3];
rz(-0.85060564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7443202) q[2];
sx q[2];
rz(-2.263415) q[2];
sx q[2];
rz(1.6142023) q[2];
rz(0.39198908) q[3];
sx q[3];
rz(-1.9422928) q[3];
sx q[3];
rz(-0.049662445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.024260661) q[0];
sx q[0];
rz(-1.4735104) q[0];
sx q[0];
rz(0.20183739) q[0];
rz(0.2233389) q[1];
sx q[1];
rz(-0.52450648) q[1];
sx q[1];
rz(-2.7490659) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.057923) q[0];
sx q[0];
rz(-1.4652848) q[0];
sx q[0];
rz(2.0353464) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3074218) q[2];
sx q[2];
rz(-1.6096707) q[2];
sx q[2];
rz(0.72134774) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9397247) q[1];
sx q[1];
rz(-1.0547259) q[1];
sx q[1];
rz(-0.49560541) q[1];
x q[2];
rz(-0.5397615) q[3];
sx q[3];
rz(-0.54815147) q[3];
sx q[3];
rz(2.2659311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4203804) q[2];
sx q[2];
rz(-2.1537697) q[2];
sx q[2];
rz(-0.58050275) q[2];
rz(0.37603363) q[3];
sx q[3];
rz(-2.1707363) q[3];
sx q[3];
rz(-1.6002801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9911983) q[0];
sx q[0];
rz(-1.5491485) q[0];
sx q[0];
rz(-1.3843672) q[0];
rz(1.4906384) q[1];
sx q[1];
rz(-1.2532267) q[1];
sx q[1];
rz(0.86029235) q[1];
rz(2.5096624) q[2];
sx q[2];
rz(-0.9117374) q[2];
sx q[2];
rz(0.82814817) q[2];
rz(1.5483472) q[3];
sx q[3];
rz(-1.2590564) q[3];
sx q[3];
rz(1.1127478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
