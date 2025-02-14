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
rz(0.97857082) q[0];
sx q[0];
rz(-0.11712722) q[0];
sx q[0];
rz(0.22476619) q[0];
rz(2.6449142) q[1];
sx q[1];
rz(4.7159046) q[1];
sx q[1];
rz(11.386303) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4311572) q[0];
sx q[0];
rz(-1.4340766) q[0];
sx q[0];
rz(1.4275292) q[0];
x q[1];
rz(1.5772555) q[2];
sx q[2];
rz(-2.0142372) q[2];
sx q[2];
rz(-0.214516) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2121316) q[1];
sx q[1];
rz(-0.60570215) q[1];
sx q[1];
rz(-2.6462161) q[1];
rz(-1.1389569) q[3];
sx q[3];
rz(-1.6434998) q[3];
sx q[3];
rz(1.7464701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.67285377) q[2];
sx q[2];
rz(-1.8127952) q[2];
sx q[2];
rz(1.2057745) q[2];
rz(-0.85637158) q[3];
sx q[3];
rz(-2.2090293) q[3];
sx q[3];
rz(2.2642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67742753) q[0];
sx q[0];
rz(-2.7441315) q[0];
sx q[0];
rz(-2.933266) q[0];
rz(-1.9147929) q[1];
sx q[1];
rz(-1.0700763) q[1];
sx q[1];
rz(-2.5770381) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0923742) q[0];
sx q[0];
rz(-2.9764247) q[0];
sx q[0];
rz(1.8333927) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6001098) q[2];
sx q[2];
rz(-0.80831203) q[2];
sx q[2];
rz(-1.186651) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37988389) q[1];
sx q[1];
rz(-1.392388) q[1];
sx q[1];
rz(-0.70990035) q[1];
rz(-pi) q[2];
rz(-1.0309593) q[3];
sx q[3];
rz(-1.5542277) q[3];
sx q[3];
rz(1.709721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63689256) q[2];
sx q[2];
rz(-1.8919287) q[2];
sx q[2];
rz(1.6553817) q[2];
rz(2.6460323) q[3];
sx q[3];
rz(-1.4080181) q[3];
sx q[3];
rz(1.0068309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9946852) q[0];
sx q[0];
rz(-0.72145307) q[0];
sx q[0];
rz(1.7578693) q[0];
rz(-1.2184628) q[1];
sx q[1];
rz(-1.0885025) q[1];
sx q[1];
rz(-2.7168435) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61439862) q[0];
sx q[0];
rz(-1.5081811) q[0];
sx q[0];
rz(-1.6024496) q[0];
x q[1];
rz(-2.9016657) q[2];
sx q[2];
rz(-1.8169962) q[2];
sx q[2];
rz(-2.4407516) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.16097611) q[1];
sx q[1];
rz(-2.7572639) q[1];
sx q[1];
rz(-0.20341441) q[1];
x q[2];
rz(2.5249486) q[3];
sx q[3];
rz(-1.5381406) q[3];
sx q[3];
rz(-0.045696229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72033) q[2];
sx q[2];
rz(-2.0739372) q[2];
sx q[2];
rz(0.024451582) q[2];
rz(-1.1484185) q[3];
sx q[3];
rz(-1.0992071) q[3];
sx q[3];
rz(1.0816157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20632437) q[0];
sx q[0];
rz(-0.22980389) q[0];
sx q[0];
rz(-2.9277053) q[0];
rz(-2.3230486) q[1];
sx q[1];
rz(-1.7844424) q[1];
sx q[1];
rz(-2.4549386) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.548107) q[0];
sx q[0];
rz(-1.0229551) q[0];
sx q[0];
rz(-2.0440897) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29149194) q[2];
sx q[2];
rz(-2.1003869) q[2];
sx q[2];
rz(-1.6232217) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80685341) q[1];
sx q[1];
rz(-2.377393) q[1];
sx q[1];
rz(3.0055226) q[1];
x q[2];
rz(-2.7044933) q[3];
sx q[3];
rz(-1.0954787) q[3];
sx q[3];
rz(-0.10224414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.116918) q[2];
sx q[2];
rz(-1.7322098) q[2];
sx q[2];
rz(0.68946687) q[2];
rz(1.3245448) q[3];
sx q[3];
rz(-1.2169415) q[3];
sx q[3];
rz(2.569258) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6784742) q[0];
sx q[0];
rz(-0.093570396) q[0];
sx q[0];
rz(2.1043188) q[0];
rz(0.67266881) q[1];
sx q[1];
rz(-1.9692407) q[1];
sx q[1];
rz(0.79576463) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0908431) q[0];
sx q[0];
rz(-1.2552869) q[0];
sx q[0];
rz(1.0361703) q[0];
rz(1.1188385) q[2];
sx q[2];
rz(-0.6217494) q[2];
sx q[2];
rz(-2.5202519) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2356253) q[1];
sx q[1];
rz(-0.86938953) q[1];
sx q[1];
rz(-2.7011306) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7618587) q[3];
sx q[3];
rz(-0.70443166) q[3];
sx q[3];
rz(1.3317127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.98046389) q[2];
sx q[2];
rz(-2.0420065) q[2];
sx q[2];
rz(0.83578342) q[2];
rz(2.8849844) q[3];
sx q[3];
rz(-2.0308688) q[3];
sx q[3];
rz(-0.49032828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68818727) q[0];
sx q[0];
rz(-1.1510993) q[0];
sx q[0];
rz(-1.2123464) q[0];
rz(1.0351828) q[1];
sx q[1];
rz(-1.4915219) q[1];
sx q[1];
rz(-2.129668) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31549227) q[0];
sx q[0];
rz(-1.4619215) q[0];
sx q[0];
rz(-1.3988162) q[0];
x q[1];
rz(1.5359319) q[2];
sx q[2];
rz(-0.30403411) q[2];
sx q[2];
rz(-2.3140098) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.052421817) q[1];
sx q[1];
rz(-0.60985205) q[1];
sx q[1];
rz(-2.3754041) q[1];
x q[2];
rz(-0.38326855) q[3];
sx q[3];
rz(-1.3196527) q[3];
sx q[3];
rz(2.0519837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6637491) q[2];
sx q[2];
rz(-2.9426844) q[2];
sx q[2];
rz(-2.193023) q[2];
rz(-0.47032022) q[3];
sx q[3];
rz(-2.0430653) q[3];
sx q[3];
rz(1.2547913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5941641) q[0];
sx q[0];
rz(-1.0110039) q[0];
sx q[0];
rz(-1.8773361) q[0];
rz(0.7803548) q[1];
sx q[1];
rz(-0.64096132) q[1];
sx q[1];
rz(2.9497214) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9148798) q[0];
sx q[0];
rz(-2.5396945) q[0];
sx q[0];
rz(-0.65897091) q[0];
rz(-pi) q[1];
rz(-2.3290995) q[2];
sx q[2];
rz(-2.6408962) q[2];
sx q[2];
rz(-2.814584) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3873008) q[1];
sx q[1];
rz(-2.6151492) q[1];
sx q[1];
rz(0.5443046) q[1];
rz(-1.0818008) q[3];
sx q[3];
rz(-2.4823175) q[3];
sx q[3];
rz(-0.47283543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.65036217) q[2];
sx q[2];
rz(-2.0872842) q[2];
sx q[2];
rz(0.80797705) q[2];
rz(-2.2982277) q[3];
sx q[3];
rz(-1.1542412) q[3];
sx q[3];
rz(-1.5850867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4845487) q[0];
sx q[0];
rz(-3.1175685) q[0];
sx q[0];
rz(2.3727544) q[0];
rz(-2.8688042) q[1];
sx q[1];
rz(-1.6285248) q[1];
sx q[1];
rz(-1.3974894) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13279877) q[0];
sx q[0];
rz(-1.9340314) q[0];
sx q[0];
rz(1.2433267) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0553841) q[2];
sx q[2];
rz(-1.2920818) q[2];
sx q[2];
rz(0.14794042) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7556909) q[1];
sx q[1];
rz(-3.0100757) q[1];
sx q[1];
rz(1.0908068) q[1];
rz(2.6121749) q[3];
sx q[3];
rz(-2.4000958) q[3];
sx q[3];
rz(1.3138258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2046795) q[2];
sx q[2];
rz(-2.0759373) q[2];
sx q[2];
rz(1.5020471) q[2];
rz(-1.8223193) q[3];
sx q[3];
rz(-2.2917512) q[3];
sx q[3];
rz(-1.5892861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(2.4424292) q[0];
sx q[0];
rz(-2.2971239) q[0];
sx q[0];
rz(2.1635639) q[0];
rz(0.49081048) q[1];
sx q[1];
rz(-1.4421137) q[1];
sx q[1];
rz(-1.4960272) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1745047) q[0];
sx q[0];
rz(-1.3347111) q[0];
sx q[0];
rz(1.4437136) q[0];
rz(-pi) q[1];
rz(0.37668682) q[2];
sx q[2];
rz(-0.44038195) q[2];
sx q[2];
rz(-0.96681606) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2435808) q[1];
sx q[1];
rz(-1.812938) q[1];
sx q[1];
rz(0.33278521) q[1];
x q[2];
rz(3.0127527) q[3];
sx q[3];
rz(-2.6059756) q[3];
sx q[3];
rz(-0.72112668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8533123) q[2];
sx q[2];
rz(-0.86709443) q[2];
sx q[2];
rz(1.8193998) q[2];
rz(-0.36618048) q[3];
sx q[3];
rz(-1.4270695) q[3];
sx q[3];
rz(-2.0100994) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2243097) q[0];
sx q[0];
rz(-1.5320822) q[0];
sx q[0];
rz(0.073534615) q[0];
rz(-1.3167943) q[1];
sx q[1];
rz(-1.2578745) q[1];
sx q[1];
rz(-1.8427461) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.038485) q[0];
sx q[0];
rz(-1.847205) q[0];
sx q[0];
rz(-1.3871219) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8180979) q[2];
sx q[2];
rz(-2.6404233) q[2];
sx q[2];
rz(2.4275818) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.21231743) q[1];
sx q[1];
rz(-1.3354509) q[1];
sx q[1];
rz(1.9566243) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73963005) q[3];
sx q[3];
rz(-1.4907903) q[3];
sx q[3];
rz(-2.7440939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9650044) q[2];
sx q[2];
rz(-2.1829288) q[2];
sx q[2];
rz(-0.57903543) q[2];
rz(-1.6252919) q[3];
sx q[3];
rz(-2.5987891) q[3];
sx q[3];
rz(1.6453086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2753006) q[0];
sx q[0];
rz(-2.0949114) q[0];
sx q[0];
rz(3.0770009) q[0];
rz(-2.2304089) q[1];
sx q[1];
rz(-1.8235527) q[1];
sx q[1];
rz(1.8240737) q[1];
rz(1.5673135) q[2];
sx q[2];
rz(-1.9949253) q[2];
sx q[2];
rz(-0.16535769) q[2];
rz(-1.961962) q[3];
sx q[3];
rz(-0.65476553) q[3];
sx q[3];
rz(-2.7162566) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
