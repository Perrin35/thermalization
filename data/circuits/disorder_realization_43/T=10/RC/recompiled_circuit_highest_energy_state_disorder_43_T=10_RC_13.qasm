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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49591767) q[0];
sx q[0];
rz(-1.1326651) q[0];
sx q[0];
rz(-2.8854831) q[0];
rz(-pi) q[1];
rz(-2.5426504) q[2];
sx q[2];
rz(-1.1556731) q[2];
sx q[2];
rz(0.51905635) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4551022) q[1];
sx q[1];
rz(-0.89340607) q[1];
sx q[1];
rz(1.6325006) q[1];
rz(1.41134) q[3];
sx q[3];
rz(-2.660369) q[3];
sx q[3];
rz(-3.0414097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.072210463) q[2];
sx q[2];
rz(-2.78077) q[2];
sx q[2];
rz(-1.3145831) q[2];
rz(0.9032816) q[3];
sx q[3];
rz(-1.3638146) q[3];
sx q[3];
rz(0.34590736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2717993) q[0];
sx q[0];
rz(-1.485774) q[0];
sx q[0];
rz(2.2858009) q[0];
rz(0.77404147) q[1];
sx q[1];
rz(-1.1639405) q[1];
sx q[1];
rz(2.3537297) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6928394) q[0];
sx q[0];
rz(-1.3200545) q[0];
sx q[0];
rz(0.022912774) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6131975) q[2];
sx q[2];
rz(-1.9147434) q[2];
sx q[2];
rz(2.6478772) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98709244) q[1];
sx q[1];
rz(-2.5804511) q[1];
sx q[1];
rz(-0.89548703) q[1];
rz(-pi) q[2];
x q[2];
rz(2.47452) q[3];
sx q[3];
rz(-1.4049585) q[3];
sx q[3];
rz(1.1950243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9539255) q[0];
sx q[0];
rz(-1.6463771) q[0];
sx q[0];
rz(0.060610108) q[0];
rz(2.4568779) q[1];
sx q[1];
rz(-1.741332) q[1];
sx q[1];
rz(-0.33904591) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1364839) q[0];
sx q[0];
rz(-1.5915131) q[0];
sx q[0];
rz(2.5584975) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9525063) q[2];
sx q[2];
rz(-1.4497309) q[2];
sx q[2];
rz(2.6775286) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.10968929) q[1];
sx q[1];
rz(-1.2983783) q[1];
sx q[1];
rz(-1.2485571) q[1];
rz(-pi) q[2];
rz(-0.77683461) q[3];
sx q[3];
rz(-1.562444) q[3];
sx q[3];
rz(1.9296822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.48245779) q[2];
sx q[2];
rz(-1.8637916) q[2];
sx q[2];
rz(-0.47046146) q[2];
rz(-0.86236924) q[3];
sx q[3];
rz(-0.44895288) q[3];
sx q[3];
rz(-1.580015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7793133) q[0];
sx q[0];
rz(-1.1184432) q[0];
sx q[0];
rz(-0.40801868) q[0];
rz(1.7904003) q[1];
sx q[1];
rz(-1.2974757) q[1];
sx q[1];
rz(1.8203576) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7119638) q[0];
sx q[0];
rz(-0.59998673) q[0];
sx q[0];
rz(-0.55498755) q[0];
x q[1];
rz(3.0555326) q[2];
sx q[2];
rz(-2.2241908) q[2];
sx q[2];
rz(-1.6376405) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.27038747) q[1];
sx q[1];
rz(-3.0000038) q[1];
sx q[1];
rz(0.54091994) q[1];
rz(-0.77797555) q[3];
sx q[3];
rz(-2.5730592) q[3];
sx q[3];
rz(-1.83873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2230175) q[2];
sx q[2];
rz(-1.431798) q[2];
sx q[2];
rz(-2.4187386) q[2];
rz(2.2919848) q[3];
sx q[3];
rz(-1.1725715) q[3];
sx q[3];
rz(-2.1037219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-1.9549114) q[0];
sx q[0];
rz(-0.51690042) q[0];
sx q[0];
rz(-2.5816259) q[0];
rz(-2.9365183) q[1];
sx q[1];
rz(-1.2566902) q[1];
sx q[1];
rz(2.8505039) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9331074) q[0];
sx q[0];
rz(-1.7400842) q[0];
sx q[0];
rz(-0.27219682) q[0];
rz(-pi) q[1];
rz(-0.79925691) q[2];
sx q[2];
rz(-1.6148668) q[2];
sx q[2];
rz(-1.26233) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2771137) q[1];
sx q[1];
rz(-0.8833589) q[1];
sx q[1];
rz(2.8836125) q[1];
rz(0.17125968) q[3];
sx q[3];
rz(-0.75920877) q[3];
sx q[3];
rz(-1.1472792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2228955) q[2];
sx q[2];
rz(-1.4783858) q[2];
sx q[2];
rz(-2.843294) q[2];
rz(-1.9722624) q[3];
sx q[3];
rz(-0.9669286) q[3];
sx q[3];
rz(1.6303308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8339612) q[0];
sx q[0];
rz(-1.6809604) q[0];
sx q[0];
rz(0.6024012) q[0];
rz(2.7700453) q[1];
sx q[1];
rz(-1.5448152) q[1];
sx q[1];
rz(1.0317624) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69799549) q[0];
sx q[0];
rz(-1.5549608) q[0];
sx q[0];
rz(-1.3965194) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6844031) q[2];
sx q[2];
rz(-2.6962099) q[2];
sx q[2];
rz(-0.066748652) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.18373016) q[1];
sx q[1];
rz(-2.0477844) q[1];
sx q[1];
rz(-1.7836003) q[1];
rz(-pi) q[2];
rz(2.2183275) q[3];
sx q[3];
rz(-2.0482691) q[3];
sx q[3];
rz(1.6357733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3370257) q[2];
sx q[2];
rz(-1.1013384) q[2];
sx q[2];
rz(1.6229013) q[2];
rz(0.8484146) q[3];
sx q[3];
rz(-2.266326) q[3];
sx q[3];
rz(-3.1019822) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0254211) q[0];
sx q[0];
rz(-0.051055901) q[0];
sx q[0];
rz(-2.1436932) q[0];
rz(-2.8669224) q[1];
sx q[1];
rz(-0.88013595) q[1];
sx q[1];
rz(-1.3708699) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3034015) q[0];
sx q[0];
rz(-2.5104021) q[0];
sx q[0];
rz(3.1188008) q[0];
rz(2.9878658) q[2];
sx q[2];
rz(-1.9228553) q[2];
sx q[2];
rz(1.597126) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9129907) q[1];
sx q[1];
rz(-0.96768846) q[1];
sx q[1];
rz(-1.2772389) q[1];
x q[2];
rz(-0.82084772) q[3];
sx q[3];
rz(-2.2162262) q[3];
sx q[3];
rz(0.57725932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.85713282) q[2];
sx q[2];
rz(-2.1993115) q[2];
sx q[2];
rz(-0.029646309) q[2];
rz(0.61521411) q[3];
sx q[3];
rz(-1.3968202) q[3];
sx q[3];
rz(2.7992547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0520332) q[0];
sx q[0];
rz(-2.679306) q[0];
sx q[0];
rz(2.3585368) q[0];
rz(2.0217333) q[1];
sx q[1];
rz(-0.66711396) q[1];
sx q[1];
rz(-1.3979744) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92358855) q[0];
sx q[0];
rz(-2.5177885) q[0];
sx q[0];
rz(0.82360928) q[0];
rz(-pi) q[1];
rz(0.2335642) q[2];
sx q[2];
rz(-1.2937045) q[2];
sx q[2];
rz(-1.8943) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.74136342) q[1];
sx q[1];
rz(-1.5619763) q[1];
sx q[1];
rz(0.90290248) q[1];
x q[2];
rz(0.78432958) q[3];
sx q[3];
rz(-1.722986) q[3];
sx q[3];
rz(-2.1209716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8457501) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(-1.8208549) q[2];
rz(2.1218421) q[3];
sx q[3];
rz(-1.7176065) q[3];
sx q[3];
rz(-1.7216871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19972292) q[0];
sx q[0];
rz(-2.0482735) q[0];
sx q[0];
rz(2.1871908) q[0];
rz(1.9873387) q[1];
sx q[1];
rz(-2.4955165) q[1];
sx q[1];
rz(-2.0571713) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99710871) q[0];
sx q[0];
rz(-0.86291828) q[0];
sx q[0];
rz(-3.0860712) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25523369) q[2];
sx q[2];
rz(-1.298578) q[2];
sx q[2];
rz(-1.9013426) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.22218765) q[1];
sx q[1];
rz(-2.0424941) q[1];
sx q[1];
rz(3.0604048) q[1];
rz(2.9436005) q[3];
sx q[3];
rz(-2.3831203) q[3];
sx q[3];
rz(0.93659479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.083100975) q[2];
sx q[2];
rz(-0.19187555) q[2];
sx q[2];
rz(0.43711883) q[2];
rz(1.7982177) q[3];
sx q[3];
rz(-1.790204) q[3];
sx q[3];
rz(0.04537151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0992391) q[0];
sx q[0];
rz(-2.1411242) q[0];
sx q[0];
rz(-1.4105256) q[0];
rz(0.22077416) q[1];
sx q[1];
rz(-0.53123728) q[1];
sx q[1];
rz(0.9332307) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8932874) q[0];
sx q[0];
rz(-2.3961005) q[0];
sx q[0];
rz(3.1214691) q[0];
x q[1];
rz(1.2104697) q[2];
sx q[2];
rz(-2.1585352) q[2];
sx q[2];
rz(-0.34229842) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1677959) q[1];
sx q[1];
rz(-0.78404155) q[1];
sx q[1];
rz(1.7586437) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0296928) q[3];
sx q[3];
rz(-0.94116941) q[3];
sx q[3];
rz(0.96523413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.71020469) q[2];
sx q[2];
rz(-0.34222558) q[2];
sx q[2];
rz(2.5101275) q[2];
rz(0.7190052) q[3];
sx q[3];
rz(-1.7944261) q[3];
sx q[3];
rz(-2.9001111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5030293) q[0];
sx q[0];
rz(-1.4894435) q[0];
sx q[0];
rz(1.3110934) q[0];
rz(-0.44166625) q[1];
sx q[1];
rz(-1.3413981) q[1];
sx q[1];
rz(-1.9955019) q[1];
rz(1.6497816) q[2];
sx q[2];
rz(-0.38182237) q[2];
sx q[2];
rz(0.67425722) q[2];
rz(0.56959116) q[3];
sx q[3];
rz(-1.4841595) q[3];
sx q[3];
rz(-1.0356173) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
