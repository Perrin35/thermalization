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
rz(-2.4515406) q[0];
sx q[0];
rz(-2.2909988) q[0];
sx q[0];
rz(2.9414862) q[0];
rz(2.9432358) q[1];
sx q[1];
rz(-2.618572) q[1];
sx q[1];
rz(-1.3318292) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057077335) q[0];
sx q[0];
rz(-2.6383218) q[0];
sx q[0];
rz(2.0664735) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66352377) q[2];
sx q[2];
rz(-0.71394701) q[2];
sx q[2];
rz(2.6235142) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7847474) q[1];
sx q[1];
rz(-2.4618399) q[1];
sx q[1];
rz(-0.076514449) q[1];
x q[2];
rz(0.08272127) q[3];
sx q[3];
rz(-1.0961856) q[3];
sx q[3];
rz(0.079291346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0693822) q[2];
sx q[2];
rz(-0.36082265) q[2];
sx q[2];
rz(-1.8270095) q[2];
rz(-0.9032816) q[3];
sx q[3];
rz(-1.7777781) q[3];
sx q[3];
rz(-2.7956853) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2717993) q[0];
sx q[0];
rz(-1.485774) q[0];
sx q[0];
rz(0.85579175) q[0];
rz(-0.77404147) q[1];
sx q[1];
rz(-1.1639405) q[1];
sx q[1];
rz(0.78786293) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0252357) q[0];
sx q[0];
rz(-1.5929925) q[0];
sx q[0];
rz(1.8216013) q[0];
rz(1.5283952) q[2];
sx q[2];
rz(-1.9147434) q[2];
sx q[2];
rz(0.49371546) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.012177906) q[1];
sx q[1];
rz(-1.2316634) q[1];
sx q[1];
rz(2.0268834) q[1];
rz(-1.7807021) q[3];
sx q[3];
rz(-2.2271101) q[3];
sx q[3];
rz(-0.24649749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0756125) q[2];
sx q[2];
rz(-1.2300666) q[2];
sx q[2];
rz(-0.39609972) q[2];
rz(1.5710477) q[3];
sx q[3];
rz(-1.6377662) q[3];
sx q[3];
rz(-1.3554696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9539255) q[0];
sx q[0];
rz(-1.6463771) q[0];
sx q[0];
rz(0.060610108) q[0];
rz(0.68471471) q[1];
sx q[1];
rz(-1.4002607) q[1];
sx q[1];
rz(-0.33904591) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6936172) q[0];
sx q[0];
rz(-0.98784271) q[0];
sx q[0];
rz(1.5459803) q[0];
x q[1];
rz(-2.5671447) q[2];
sx q[2];
rz(-2.9174605) q[2];
sx q[2];
rz(1.4719065) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3587957) q[1];
sx q[1];
rz(-2.7227245) q[1];
sx q[1];
rz(2.2936506) q[1];
rz(-0.77683461) q[3];
sx q[3];
rz(-1.5791487) q[3];
sx q[3];
rz(1.2119105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.48245779) q[2];
sx q[2];
rz(-1.277801) q[2];
sx q[2];
rz(0.47046146) q[2];
rz(0.86236924) q[3];
sx q[3];
rz(-2.6926398) q[3];
sx q[3];
rz(-1.580015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.7793133) q[0];
sx q[0];
rz(-2.0231495) q[0];
sx q[0];
rz(-0.40801868) q[0];
rz(-1.7904003) q[1];
sx q[1];
rz(-1.2974757) q[1];
sx q[1];
rz(1.321235) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7119638) q[0];
sx q[0];
rz(-0.59998673) q[0];
sx q[0];
rz(-2.5866051) q[0];
rz(-pi) q[1];
rz(-0.91561134) q[2];
sx q[2];
rz(-1.639099) q[2];
sx q[2];
rz(-0.11924041) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3776748) q[1];
sx q[1];
rz(-1.6435247) q[1];
sx q[1];
rz(-3.0200028) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9923237) q[3];
sx q[3];
rz(-1.9644004) q[3];
sx q[3];
rz(0.43962653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2230175) q[2];
sx q[2];
rz(-1.7097946) q[2];
sx q[2];
rz(-2.4187386) q[2];
rz(2.2919848) q[3];
sx q[3];
rz(-1.1725715) q[3];
sx q[3];
rz(1.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.9549114) q[0];
sx q[0];
rz(-2.6246922) q[0];
sx q[0];
rz(2.5816259) q[0];
rz(2.9365183) q[1];
sx q[1];
rz(-1.2566902) q[1];
sx q[1];
rz(-2.8505039) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2084853) q[0];
sx q[0];
rz(-1.4015084) q[0];
sx q[0];
rz(0.27219682) q[0];
rz(3.0801513) q[2];
sx q[2];
rz(-2.3413918) q[2];
sx q[2];
rz(0.35129181) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4704125) q[1];
sx q[1];
rz(-2.414783) q[1];
sx q[1];
rz(-1.2695168) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.970333) q[3];
sx q[3];
rz(-2.3823839) q[3];
sx q[3];
rz(-1.9943135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2228955) q[2];
sx q[2];
rz(-1.6632068) q[2];
sx q[2];
rz(-0.29829868) q[2];
rz(1.9722624) q[3];
sx q[3];
rz(-2.1746641) q[3];
sx q[3];
rz(-1.5112618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30763141) q[0];
sx q[0];
rz(-1.4606322) q[0];
sx q[0];
rz(-2.5391915) q[0];
rz(2.7700453) q[1];
sx q[1];
rz(-1.5967775) q[1];
sx q[1];
rz(-1.0317624) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69799549) q[0];
sx q[0];
rz(-1.5866318) q[0];
sx q[0];
rz(1.7450733) q[0];
rz(-pi) q[1];
rz(2.736892) q[2];
sx q[2];
rz(-1.3794624) q[2];
sx q[2];
rz(-1.2197242) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9578625) q[1];
sx q[1];
rz(-1.0938083) q[1];
sx q[1];
rz(-1.3579923) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5661181) q[3];
sx q[3];
rz(-2.1361975) q[3];
sx q[3];
rz(2.7421212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3370257) q[2];
sx q[2];
rz(-1.1013384) q[2];
sx q[2];
rz(1.6229013) q[2];
rz(2.2931781) q[3];
sx q[3];
rz(-0.87526667) q[3];
sx q[3];
rz(0.039610473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0254211) q[0];
sx q[0];
rz(-3.0905368) q[0];
sx q[0];
rz(2.1436932) q[0];
rz(-2.8669224) q[1];
sx q[1];
rz(-0.88013595) q[1];
sx q[1];
rz(-1.3708699) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8664198) q[0];
sx q[0];
rz(-0.93979561) q[0];
sx q[0];
rz(-1.5541398) q[0];
rz(-pi) q[1];
rz(-2.9878658) q[2];
sx q[2];
rz(-1.9228553) q[2];
sx q[2];
rz(-1.597126) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8805367) q[1];
sx q[1];
rz(-2.4789101) q[1];
sx q[1];
rz(2.7438394) q[1];
rz(-pi) q[2];
rz(2.4059345) q[3];
sx q[3];
rz(-2.1950985) q[3];
sx q[3];
rz(-1.566938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.85713282) q[2];
sx q[2];
rz(-0.94228116) q[2];
sx q[2];
rz(-3.1119463) q[2];
rz(-2.5263785) q[3];
sx q[3];
rz(-1.3968202) q[3];
sx q[3];
rz(2.7992547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089559473) q[0];
sx q[0];
rz(-2.679306) q[0];
sx q[0];
rz(0.7830559) q[0];
rz(-1.1198593) q[1];
sx q[1];
rz(-2.4744787) q[1];
sx q[1];
rz(-1.7436183) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0692978) q[0];
sx q[0];
rz(-2.0136478) q[0];
sx q[0];
rz(-2.686727) q[0];
x q[1];
rz(2.2538848) q[2];
sx q[2];
rz(-2.7811433) q[2];
sx q[2];
rz(0.53149022) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3009792) q[1];
sx q[1];
rz(-2.4736495) q[1];
sx q[1];
rz(1.5850369) q[1];
x q[2];
rz(-2.9277738) q[3];
sx q[3];
rz(-0.79584661) q[3];
sx q[3];
rz(0.70094943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.29584259) q[2];
sx q[2];
rz(-2.5551899) q[2];
sx q[2];
rz(-1.8208549) q[2];
rz(1.0197506) q[3];
sx q[3];
rz(-1.7176065) q[3];
sx q[3];
rz(-1.4199055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9418697) q[0];
sx q[0];
rz(-1.0933192) q[0];
sx q[0];
rz(2.1871908) q[0];
rz(1.154254) q[1];
sx q[1];
rz(-0.64607611) q[1];
sx q[1];
rz(-2.0571713) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1444839) q[0];
sx q[0];
rz(-0.86291828) q[0];
sx q[0];
rz(-3.0860712) q[0];
x q[1];
rz(1.8516638) q[2];
sx q[2];
rz(-1.8164338) q[2];
sx q[2];
rz(-0.40058595) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0965914) q[1];
sx q[1];
rz(-0.47811478) q[1];
sx q[1];
rz(1.4131312) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7550657) q[3];
sx q[3];
rz(-2.3109155) q[3];
sx q[3];
rz(-0.6669464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.083100975) q[2];
sx q[2];
rz(-0.19187555) q[2];
sx q[2];
rz(0.43711883) q[2];
rz(-1.7982177) q[3];
sx q[3];
rz(-1.3513887) q[3];
sx q[3];
rz(0.04537151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0992391) q[0];
sx q[0];
rz(-1.0004685) q[0];
sx q[0];
rz(1.7310671) q[0];
rz(-0.22077416) q[1];
sx q[1];
rz(-0.53123728) q[1];
sx q[1];
rz(-0.9332307) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2756903) q[0];
sx q[0];
rz(-0.82549107) q[0];
sx q[0];
rz(1.5893713) q[0];
x q[1];
rz(-2.5228517) q[2];
sx q[2];
rz(-1.2729984) q[2];
sx q[2];
rz(2.1190475) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9737967) q[1];
sx q[1];
rz(-0.78404155) q[1];
sx q[1];
rz(1.7586437) q[1];
x q[2];
rz(0.93818061) q[3];
sx q[3];
rz(-1.6611735) q[3];
sx q[3];
rz(2.602102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.431388) q[2];
sx q[2];
rz(-0.34222558) q[2];
sx q[2];
rz(0.63146511) q[2];
rz(-2.4225875) q[3];
sx q[3];
rz(-1.7944261) q[3];
sx q[3];
rz(0.24148153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5030293) q[0];
sx q[0];
rz(-1.4894435) q[0];
sx q[0];
rz(1.3110934) q[0];
rz(2.6999264) q[1];
sx q[1];
rz(-1.3413981) q[1];
sx q[1];
rz(-1.9955019) q[1];
rz(-1.491811) q[2];
sx q[2];
rz(-0.38182237) q[2];
sx q[2];
rz(0.67425722) q[2];
rz(-0.15968237) q[3];
sx q[3];
rz(-2.5661712) q[3];
sx q[3];
rz(-2.4721091) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
