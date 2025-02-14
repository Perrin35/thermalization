OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0491068) q[0];
sx q[0];
rz(-1.8422814) q[0];
sx q[0];
rz(0.045825034) q[0];
rz(1.2996281) q[1];
sx q[1];
rz(-1.8027432) q[1];
sx q[1];
rz(1.2531228) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50025138) q[0];
sx q[0];
rz(-2.0098899) q[0];
sx q[0];
rz(2.8706949) q[0];
rz(-pi) q[1];
rz(1.5270751) q[2];
sx q[2];
rz(-1.7720024) q[2];
sx q[2];
rz(-1.218856) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.040267) q[1];
sx q[1];
rz(-0.76124746) q[1];
sx q[1];
rz(2.9684307) q[1];
rz(0.0039611749) q[3];
sx q[3];
rz(-1.4433493) q[3];
sx q[3];
rz(-2.733903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32156285) q[2];
sx q[2];
rz(-1.4048615) q[2];
sx q[2];
rz(-0.31537867) q[2];
rz(-3.0553715) q[3];
sx q[3];
rz(-0.092970522) q[3];
sx q[3];
rz(-0.84403795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.740199) q[0];
sx q[0];
rz(-1.711015) q[0];
sx q[0];
rz(-2.8061818) q[0];
rz(-0.3849349) q[1];
sx q[1];
rz(-0.79610577) q[1];
sx q[1];
rz(-1.8523432) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1803363) q[0];
sx q[0];
rz(-1.633005) q[0];
sx q[0];
rz(0.38952413) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.697549) q[2];
sx q[2];
rz(-0.64536649) q[2];
sx q[2];
rz(0.54189787) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1114167) q[1];
sx q[1];
rz(-1.9523638) q[1];
sx q[1];
rz(-1.2658582) q[1];
rz(-1.3277169) q[3];
sx q[3];
rz(-1.7020924) q[3];
sx q[3];
rz(2.8158902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.59911597) q[2];
sx q[2];
rz(-1.3834388) q[2];
sx q[2];
rz(1.4625589) q[2];
rz(-0.8820495) q[3];
sx q[3];
rz(-0.65989152) q[3];
sx q[3];
rz(1.5006458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3682692) q[0];
sx q[0];
rz(-0.48650807) q[0];
sx q[0];
rz(0.64273709) q[0];
rz(0.87359387) q[1];
sx q[1];
rz(-1.3106376) q[1];
sx q[1];
rz(2.1547623) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5499134) q[0];
sx q[0];
rz(-2.4351623) q[0];
sx q[0];
rz(-1.3384992) q[0];
rz(-2.0552733) q[2];
sx q[2];
rz(-1.6256126) q[2];
sx q[2];
rz(-1.9475186) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.75007406) q[1];
sx q[1];
rz(-2.1100535) q[1];
sx q[1];
rz(0.92309322) q[1];
x q[2];
rz(1.6968459) q[3];
sx q[3];
rz(-2.7641962) q[3];
sx q[3];
rz(-2.6871329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.70283908) q[2];
sx q[2];
rz(-0.025391014) q[2];
sx q[2];
rz(-1.0566443) q[2];
rz(1.7993401) q[3];
sx q[3];
rz(-1.4547576) q[3];
sx q[3];
rz(-0.87856436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.037828) q[0];
sx q[0];
rz(-0.28452888) q[0];
sx q[0];
rz(-1.9985265) q[0];
rz(0.68483886) q[1];
sx q[1];
rz(-1.7999444) q[1];
sx q[1];
rz(0.24982223) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046127675) q[0];
sx q[0];
rz(-0.85554142) q[0];
sx q[0];
rz(-0.94761316) q[0];
rz(-pi) q[1];
rz(-0.5742214) q[2];
sx q[2];
rz(-1.5674233) q[2];
sx q[2];
rz(-2.8115438) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4433684) q[1];
sx q[1];
rz(-2.0495546) q[1];
sx q[1];
rz(2.4633292) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3288199) q[3];
sx q[3];
rz(-1.5723088) q[3];
sx q[3];
rz(-0.63336271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.41986856) q[2];
sx q[2];
rz(-1.3519752) q[2];
sx q[2];
rz(1.8458337) q[2];
rz(-1.7660247) q[3];
sx q[3];
rz(-1.7121168) q[3];
sx q[3];
rz(-3.0425369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38882035) q[0];
sx q[0];
rz(-1.7701912) q[0];
sx q[0];
rz(-0.8465299) q[0];
rz(1.9580152) q[1];
sx q[1];
rz(-1.9639683) q[1];
sx q[1];
rz(1.2394946) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.469231) q[0];
sx q[0];
rz(-1.8235908) q[0];
sx q[0];
rz(-0.12427788) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7527649) q[2];
sx q[2];
rz(-1.2335868) q[2];
sx q[2];
rz(-2.6621001) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.073055) q[1];
sx q[1];
rz(-3.0843711) q[1];
sx q[1];
rz(-1.7806609) q[1];
x q[2];
rz(-2.1679513) q[3];
sx q[3];
rz(-1.5212562) q[3];
sx q[3];
rz(1.7907536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5014629) q[2];
sx q[2];
rz(-1.632246) q[2];
sx q[2];
rz(2.787369) q[2];
rz(-0.7192449) q[3];
sx q[3];
rz(-2.5651599) q[3];
sx q[3];
rz(2.332212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7648014) q[0];
sx q[0];
rz(-2.9053595) q[0];
sx q[0];
rz(0.34673044) q[0];
rz(-2.4193343) q[1];
sx q[1];
rz(-1.6112695) q[1];
sx q[1];
rz(-0.17682704) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2472025) q[0];
sx q[0];
rz(-2.6567334) q[0];
sx q[0];
rz(-1.8762183) q[0];
x q[1];
rz(1.7406101) q[2];
sx q[2];
rz(-0.50504518) q[2];
sx q[2];
rz(1.1894047) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36784259) q[1];
sx q[1];
rz(-0.23077877) q[1];
sx q[1];
rz(-2.6683067) q[1];
rz(-1.4106287) q[3];
sx q[3];
rz(-0.56863848) q[3];
sx q[3];
rz(1.7344432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5727545) q[2];
sx q[2];
rz(-0.73358959) q[2];
sx q[2];
rz(2.0983762) q[2];
rz(-2.9853232) q[3];
sx q[3];
rz(-1.0633435) q[3];
sx q[3];
rz(3.0472896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25065502) q[0];
sx q[0];
rz(-1.0303048) q[0];
sx q[0];
rz(1.2283196) q[0];
rz(0.43371513) q[1];
sx q[1];
rz(-2.3408196) q[1];
sx q[1];
rz(1.0450276) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50253326) q[0];
sx q[0];
rz(-0.21083388) q[0];
sx q[0];
rz(1.2913713) q[0];
rz(-pi) q[1];
rz(-2.074998) q[2];
sx q[2];
rz(-0.81461755) q[2];
sx q[2];
rz(2.2058918) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.24416395) q[1];
sx q[1];
rz(-0.62627316) q[1];
sx q[1];
rz(-1.3757964) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77680334) q[3];
sx q[3];
rz(-1.9597199) q[3];
sx q[3];
rz(-0.17059205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.23994437) q[2];
sx q[2];
rz(-2.1176691) q[2];
sx q[2];
rz(-2.955692) q[2];
rz(1.9013885) q[3];
sx q[3];
rz(-1.5408206) q[3];
sx q[3];
rz(2.127229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1321201) q[0];
sx q[0];
rz(-1.8931696) q[0];
sx q[0];
rz(-0.15952071) q[0];
rz(-0.80687833) q[1];
sx q[1];
rz(-1.2346377) q[1];
sx q[1];
rz(0.0085011403) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5497564) q[0];
sx q[0];
rz(-1.4167794) q[0];
sx q[0];
rz(-1.9817673) q[0];
rz(-pi) q[1];
rz(0.79719724) q[2];
sx q[2];
rz(-1.6455499) q[2];
sx q[2];
rz(-2.9644186) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9125625) q[1];
sx q[1];
rz(-2.8797315) q[1];
sx q[1];
rz(0.53066855) q[1];
rz(-1.0866585) q[3];
sx q[3];
rz(-1.7268506) q[3];
sx q[3];
rz(-0.54359667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27935091) q[2];
sx q[2];
rz(-1.778435) q[2];
sx q[2];
rz(-0.19374338) q[2];
rz(-1.6019609) q[3];
sx q[3];
rz(-1.1774457) q[3];
sx q[3];
rz(1.1295454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9047456) q[0];
sx q[0];
rz(-1.3847677) q[0];
sx q[0];
rz(-2.4054476) q[0];
rz(-2.2537117) q[1];
sx q[1];
rz(-0.45279756) q[1];
sx q[1];
rz(-0.050051659) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0808635) q[0];
sx q[0];
rz(-1.2300778) q[0];
sx q[0];
rz(-1.3843892) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4718945) q[2];
sx q[2];
rz(-1.8840232) q[2];
sx q[2];
rz(-3.1196496) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60683317) q[1];
sx q[1];
rz(-1.8799361) q[1];
sx q[1];
rz(1.3571795) q[1];
x q[2];
rz(-0.83733674) q[3];
sx q[3];
rz(-2.3792107) q[3];
sx q[3];
rz(1.9782982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8686409) q[2];
sx q[2];
rz(-1.4401888) q[2];
sx q[2];
rz(-3.1325373) q[2];
rz(-0.35033369) q[3];
sx q[3];
rz(-0.30600268) q[3];
sx q[3];
rz(2.8370324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7250799) q[0];
sx q[0];
rz(-1.059499) q[0];
sx q[0];
rz(0.9730202) q[0];
rz(-1.579772) q[1];
sx q[1];
rz(-0.90619722) q[1];
sx q[1];
rz(0.80950338) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3448787) q[0];
sx q[0];
rz(-2.2911706) q[0];
sx q[0];
rz(1.0864837) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8000748) q[2];
sx q[2];
rz(-1.7588758) q[2];
sx q[2];
rz(-3.0517202) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1228408) q[1];
sx q[1];
rz(-2.0301308) q[1];
sx q[1];
rz(0.18045119) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91416759) q[3];
sx q[3];
rz(-2.4501213) q[3];
sx q[3];
rz(-1.4007614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1805264) q[2];
sx q[2];
rz(-0.82225353) q[2];
sx q[2];
rz(-2.6194438) q[2];
rz(-0.3565878) q[3];
sx q[3];
rz(-2.5643189) q[3];
sx q[3];
rz(-3.0780415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4840354) q[0];
sx q[0];
rz(-0.35150305) q[0];
sx q[0];
rz(1.2765314) q[0];
rz(-2.7784078) q[1];
sx q[1];
rz(-0.52774944) q[1];
sx q[1];
rz(-1.48762) q[1];
rz(2.1073916) q[2];
sx q[2];
rz(-0.86884872) q[2];
sx q[2];
rz(1.6128513) q[2];
rz(0.47388062) q[3];
sx q[3];
rz(-2.5070179) q[3];
sx q[3];
rz(-1.0567155) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
