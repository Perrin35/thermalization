OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.428838) q[0];
sx q[0];
rz(-1.5313671) q[0];
sx q[0];
rz(2.03696) q[0];
rz(1.334335) q[1];
sx q[1];
rz(-2.4185138) q[1];
sx q[1];
rz(-1.2637957) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6440463) q[0];
sx q[0];
rz(-0.9263557) q[0];
sx q[0];
rz(0.30797663) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42277067) q[2];
sx q[2];
rz(-2.1029933) q[2];
sx q[2];
rz(-1.339004) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7967148) q[1];
sx q[1];
rz(-1.895322) q[1];
sx q[1];
rz(3.0561563) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0588082) q[3];
sx q[3];
rz(-0.48824874) q[3];
sx q[3];
rz(1.4940718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9870712) q[2];
sx q[2];
rz(-1.0769341) q[2];
sx q[2];
rz(-1.3686352) q[2];
rz(-0.23855071) q[3];
sx q[3];
rz(-2.7261901) q[3];
sx q[3];
rz(0.11428782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.7411165) q[0];
sx q[0];
rz(-1.1392765) q[0];
sx q[0];
rz(1.9912632) q[0];
rz(-0.087609619) q[1];
sx q[1];
rz(-1.3299512) q[1];
sx q[1];
rz(-1.5709343) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8920736) q[0];
sx q[0];
rz(-0.49580916) q[0];
sx q[0];
rz(0.22479381) q[0];
rz(-pi) q[1];
rz(-2.803844) q[2];
sx q[2];
rz(-0.18924533) q[2];
sx q[2];
rz(-2.7836329) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1521235) q[1];
sx q[1];
rz(-1.996576) q[1];
sx q[1];
rz(0.87321891) q[1];
rz(-2.3493166) q[3];
sx q[3];
rz(-1.2229007) q[3];
sx q[3];
rz(2.874305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.34438434) q[2];
sx q[2];
rz(-0.71343652) q[2];
sx q[2];
rz(2.9551282) q[2];
rz(-2.3421085) q[3];
sx q[3];
rz(-1.5954285) q[3];
sx q[3];
rz(-0.98108393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28476533) q[0];
sx q[0];
rz(-1.7600049) q[0];
sx q[0];
rz(0.10936603) q[0];
rz(0.60375396) q[1];
sx q[1];
rz(-1.3394638) q[1];
sx q[1];
rz(2.107479) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4277735) q[0];
sx q[0];
rz(-1.3610916) q[0];
sx q[0];
rz(-3.0056277) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66884235) q[2];
sx q[2];
rz(-0.1136264) q[2];
sx q[2];
rz(-2.4955028) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8960524) q[1];
sx q[1];
rz(-2.8421092) q[1];
sx q[1];
rz(-0.36548945) q[1];
rz(-pi) q[2];
rz(0.14520653) q[3];
sx q[3];
rz(-2.3646486) q[3];
sx q[3];
rz(0.29175419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5321396) q[2];
sx q[2];
rz(-0.9708465) q[2];
sx q[2];
rz(-2.5965221) q[2];
rz(-2.3405781) q[3];
sx q[3];
rz(-0.89371926) q[3];
sx q[3];
rz(0.092122294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6545749) q[0];
sx q[0];
rz(-1.4734522) q[0];
sx q[0];
rz(-2.6825478) q[0];
rz(2.0252939) q[1];
sx q[1];
rz(-2.3479159) q[1];
sx q[1];
rz(2.3150516) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61454489) q[0];
sx q[0];
rz(-1.2189208) q[0];
sx q[0];
rz(2.0195168) q[0];
rz(2.3514296) q[2];
sx q[2];
rz(-3.003267) q[2];
sx q[2];
rz(-2.9776855) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3459119) q[1];
sx q[1];
rz(-1.7670146) q[1];
sx q[1];
rz(-2.1225342) q[1];
rz(0.32966726) q[3];
sx q[3];
rz(-2.7894944) q[3];
sx q[3];
rz(0.47131495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35215846) q[2];
sx q[2];
rz(-2.9106079) q[2];
sx q[2];
rz(-2.3918772) q[2];
rz(-2.0189144) q[3];
sx q[3];
rz(-1.2238945) q[3];
sx q[3];
rz(2.1813755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4753251) q[0];
sx q[0];
rz(-2.1551977) q[0];
sx q[0];
rz(3.0295897) q[0];
rz(-0.22459596) q[1];
sx q[1];
rz(-1.1677531) q[1];
sx q[1];
rz(-2.3602643) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79271475) q[0];
sx q[0];
rz(-1.9747866) q[0];
sx q[0];
rz(-3.1125665) q[0];
x q[1];
rz(-1.4766598) q[2];
sx q[2];
rz(-1.899666) q[2];
sx q[2];
rz(-0.48966792) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0722149) q[1];
sx q[1];
rz(-2.2615914) q[1];
sx q[1];
rz(-2.1363972) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.16181337) q[3];
sx q[3];
rz(-2.1573632) q[3];
sx q[3];
rz(-0.88255054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3096699) q[2];
sx q[2];
rz(-0.90193844) q[2];
sx q[2];
rz(0.93264467) q[2];
rz(-2.0617088) q[3];
sx q[3];
rz(-2.6974758) q[3];
sx q[3];
rz(-0.035695765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90750736) q[0];
sx q[0];
rz(-1.1312753) q[0];
sx q[0];
rz(1.3908516) q[0];
rz(-2.8098409) q[1];
sx q[1];
rz(-1.2650047) q[1];
sx q[1];
rz(1.5501685) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46813477) q[0];
sx q[0];
rz(-0.64366699) q[0];
sx q[0];
rz(1.629384) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6076902) q[2];
sx q[2];
rz(-1.1114745) q[2];
sx q[2];
rz(1.5609891) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0968714) q[1];
sx q[1];
rz(-1.7270178) q[1];
sx q[1];
rz(0.12466431) q[1];
x q[2];
rz(-0.36339001) q[3];
sx q[3];
rz(-1.9295606) q[3];
sx q[3];
rz(3.0409558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3127689) q[2];
sx q[2];
rz(-2.2569816) q[2];
sx q[2];
rz(1.620232) q[2];
rz(-0.96794266) q[3];
sx q[3];
rz(-1.5588375) q[3];
sx q[3];
rz(2.2255285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62436002) q[0];
sx q[0];
rz(-2.7150798) q[0];
sx q[0];
rz(0.17159167) q[0];
rz(2.102237) q[1];
sx q[1];
rz(-1.2587222) q[1];
sx q[1];
rz(-1.4412057) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2644669) q[0];
sx q[0];
rz(-1.6949589) q[0];
sx q[0];
rz(-0.32372253) q[0];
rz(-1.6827379) q[2];
sx q[2];
rz(-1.0575235) q[2];
sx q[2];
rz(-1.3789195) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6771779) q[1];
sx q[1];
rz(-2.5485793) q[1];
sx q[1];
rz(-2.2423184) q[1];
x q[2];
rz(-1.7923545) q[3];
sx q[3];
rz(-2.5668813) q[3];
sx q[3];
rz(-1.2971085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8695996) q[2];
sx q[2];
rz(-1.4078434) q[2];
sx q[2];
rz(-0.54455152) q[2];
rz(-0.49992391) q[3];
sx q[3];
rz(-1.0285503) q[3];
sx q[3];
rz(-2.669529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8498103) q[0];
sx q[0];
rz(-2.6047459) q[0];
sx q[0];
rz(-0.52919069) q[0];
rz(-2.5657907) q[1];
sx q[1];
rz(-1.878783) q[1];
sx q[1];
rz(-2.0226488) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7125583) q[0];
sx q[0];
rz(-3.0858485) q[0];
sx q[0];
rz(-0.38614614) q[0];
rz(-pi) q[1];
rz(2.697102) q[2];
sx q[2];
rz(-0.40676446) q[2];
sx q[2];
rz(-0.14061804) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8994979) q[1];
sx q[1];
rz(-2.3589954) q[1];
sx q[1];
rz(0.67983869) q[1];
rz(-pi) q[2];
rz(-2.0980074) q[3];
sx q[3];
rz(-1.7526748) q[3];
sx q[3];
rz(-1.7480237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90073663) q[2];
sx q[2];
rz(-2.6764328) q[2];
sx q[2];
rz(-0.71211234) q[2];
rz(1.5322878) q[3];
sx q[3];
rz(-0.94523793) q[3];
sx q[3];
rz(1.3473264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79494548) q[0];
sx q[0];
rz(-1.1279339) q[0];
sx q[0];
rz(2.0704863) q[0];
rz(1.935293) q[1];
sx q[1];
rz(-1.0164398) q[1];
sx q[1];
rz(1.4869022) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1893715) q[0];
sx q[0];
rz(-0.69714386) q[0];
sx q[0];
rz(2.8482262) q[0];
rz(2.5096312) q[2];
sx q[2];
rz(-0.8425396) q[2];
sx q[2];
rz(-1.8441083) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2675954) q[1];
sx q[1];
rz(-1.2800299) q[1];
sx q[1];
rz(1.3895172) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8603691) q[3];
sx q[3];
rz(-1.1373925) q[3];
sx q[3];
rz(0.12896695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.508076) q[2];
sx q[2];
rz(-0.8997007) q[2];
sx q[2];
rz(-3.0637975) q[2];
rz(0.22732321) q[3];
sx q[3];
rz(-1.1319755) q[3];
sx q[3];
rz(-2.0859065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8530497) q[0];
sx q[0];
rz(-0.74497861) q[0];
sx q[0];
rz(-2.7689834) q[0];
rz(-0.44652069) q[1];
sx q[1];
rz(-1.4563072) q[1];
sx q[1];
rz(-2.2850697) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83607446) q[0];
sx q[0];
rz(-1.627632) q[0];
sx q[0];
rz(-1.5268121) q[0];
rz(3.0356016) q[2];
sx q[2];
rz(-1.7039095) q[2];
sx q[2];
rz(-0.0078545257) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.23736542) q[1];
sx q[1];
rz(-1.6886535) q[1];
sx q[1];
rz(1.9709644) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4803745) q[3];
sx q[3];
rz(-2.2228129) q[3];
sx q[3];
rz(0.87424247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7591758) q[2];
sx q[2];
rz(-2.0903812) q[2];
sx q[2];
rz(2.5860533) q[2];
rz(-2.7555452) q[3];
sx q[3];
rz(-1.6470393) q[3];
sx q[3];
rz(-2.4773795) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1356708) q[0];
sx q[0];
rz(-1.6914524) q[0];
sx q[0];
rz(1.2941262) q[0];
rz(2.338943) q[1];
sx q[1];
rz(-1.6812656) q[1];
sx q[1];
rz(0.38801286) q[1];
rz(-2.8254208) q[2];
sx q[2];
rz(-0.36973047) q[2];
sx q[2];
rz(0.75512259) q[2];
rz(-1.4925719) q[3];
sx q[3];
rz(-0.73220069) q[3];
sx q[3];
rz(-2.2673741) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
