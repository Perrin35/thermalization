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
rz(0.053057916) q[0];
sx q[0];
rz(0.36202708) q[0];
sx q[0];
rz(10.590635) q[0];
rz(-0.8968269) q[1];
sx q[1];
rz(-1.4520175) q[1];
sx q[1];
rz(-1.7133763) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8667158) q[0];
sx q[0];
rz(-0.51565114) q[0];
sx q[0];
rz(3.0192119) q[0];
x q[1];
rz(0.27730242) q[2];
sx q[2];
rz(-2.8805974) q[2];
sx q[2];
rz(-2.9948611) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5184091) q[1];
sx q[1];
rz(-1.759638) q[1];
sx q[1];
rz(0.13407003) q[1];
rz(-pi) q[2];
rz(-2.0588097) q[3];
sx q[3];
rz(-1.8086156) q[3];
sx q[3];
rz(-2.2453515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6894138) q[2];
sx q[2];
rz(-3.1248326) q[2];
sx q[2];
rz(-0.20081271) q[2];
rz(0.14874841) q[3];
sx q[3];
rz(-0.0047618682) q[3];
sx q[3];
rz(0.31140056) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0018604) q[0];
sx q[0];
rz(-0.59356028) q[0];
sx q[0];
rz(2.1019905) q[0];
rz(0.014558583) q[1];
sx q[1];
rz(-1.2332375) q[1];
sx q[1];
rz(1.5537517) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4400875) q[0];
sx q[0];
rz(-1.0453512) q[0];
sx q[0];
rz(-2.3129025) q[0];
rz(-3.1268979) q[2];
sx q[2];
rz(-1.6449682) q[2];
sx q[2];
rz(-1.6770404) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5841516) q[1];
sx q[1];
rz(-1.5366035) q[1];
sx q[1];
rz(1.8299915) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.763299) q[3];
sx q[3];
rz(-1.8979591) q[3];
sx q[3];
rz(-3.134789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.62850922) q[2];
sx q[2];
rz(-1.6078948) q[2];
sx q[2];
rz(1.7536564) q[2];
rz(1.7657109) q[3];
sx q[3];
rz(-1.0466156) q[3];
sx q[3];
rz(-2.8468813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.3917711) q[0];
sx q[0];
rz(-0.23673683) q[0];
sx q[0];
rz(2.5340875) q[0];
rz(-1.5982184) q[1];
sx q[1];
rz(-0.18074712) q[1];
sx q[1];
rz(-0.9651331) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34244363) q[0];
sx q[0];
rz(-1.5502872) q[0];
sx q[0];
rz(3.1375225) q[0];
x q[1];
rz(-2.0320545) q[2];
sx q[2];
rz(-1.8806305) q[2];
sx q[2];
rz(-0.43333915) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7399921) q[1];
sx q[1];
rz(-1.5349746) q[1];
sx q[1];
rz(1.4261246) q[1];
x q[2];
rz(0.083250982) q[3];
sx q[3];
rz(-1.6663345) q[3];
sx q[3];
rz(2.9527309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33660108) q[2];
sx q[2];
rz(-0.6707297) q[2];
sx q[2];
rz(2.2750308) q[2];
rz(2.0349515) q[3];
sx q[3];
rz(-1.589078) q[3];
sx q[3];
rz(-1.470587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.5690145) q[0];
sx q[0];
rz(-0.67936474) q[0];
sx q[0];
rz(1.6160075) q[0];
rz(-0.010604803) q[1];
sx q[1];
rz(-0.0037825982) q[1];
sx q[1];
rz(-0.7503646) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3748659) q[0];
sx q[0];
rz(-0.9095053) q[0];
sx q[0];
rz(-1.6910716) q[0];
x q[1];
rz(2.8192873) q[2];
sx q[2];
rz(-0.63045972) q[2];
sx q[2];
rz(-0.96994627) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6746241) q[1];
sx q[1];
rz(-1.4556938) q[1];
sx q[1];
rz(2.6323167) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37639002) q[3];
sx q[3];
rz(-0.3335267) q[3];
sx q[3];
rz(0.1025478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41287199) q[2];
sx q[2];
rz(-2.0527288) q[2];
sx q[2];
rz(1.9000165) q[2];
rz(3.1359172) q[3];
sx q[3];
rz(-0.81040502) q[3];
sx q[3];
rz(-2.2976105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64549696) q[0];
sx q[0];
rz(-0.050300751) q[0];
sx q[0];
rz(0.92329931) q[0];
rz(-2.3410489) q[1];
sx q[1];
rz(-0.0034595483) q[1];
sx q[1];
rz(-0.18005767) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8220953) q[0];
sx q[0];
rz(-1.8094157) q[0];
sx q[0];
rz(-3.0988295) q[0];
rz(-0.098131432) q[2];
sx q[2];
rz(-1.4906825) q[2];
sx q[2];
rz(1.1289885) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.48985607) q[1];
sx q[1];
rz(-1.3313313) q[1];
sx q[1];
rz(-0.36291562) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11074825) q[3];
sx q[3];
rz(-1.9344653) q[3];
sx q[3];
rz(-1.2895541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.30183145) q[2];
sx q[2];
rz(-1.8415035) q[2];
sx q[2];
rz(1.4361471) q[2];
rz(1.885421) q[3];
sx q[3];
rz(-1.4830282) q[3];
sx q[3];
rz(3.0459611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.489478) q[0];
sx q[0];
rz(-2.5697932) q[0];
sx q[0];
rz(-0.44334626) q[0];
rz(1.9709142) q[1];
sx q[1];
rz(-0.0010633855) q[1];
sx q[1];
rz(0.43177691) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74878377) q[0];
sx q[0];
rz(-1.8816299) q[0];
sx q[0];
rz(0.75324599) q[0];
rz(-3.058039) q[2];
sx q[2];
rz(-1.7423034) q[2];
sx q[2];
rz(-0.87860859) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1990423) q[1];
sx q[1];
rz(-1.1909134) q[1];
sx q[1];
rz(0.6013077) q[1];
x q[2];
rz(-2.5434142) q[3];
sx q[3];
rz(-0.39387396) q[3];
sx q[3];
rz(2.1109714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40054896) q[2];
sx q[2];
rz(-2.2811175) q[2];
sx q[2];
rz(1.7575556) q[2];
rz(-0.69796973) q[3];
sx q[3];
rz(-2.3845086) q[3];
sx q[3];
rz(2.82011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31239241) q[0];
sx q[0];
rz(-1.3425403) q[0];
sx q[0];
rz(-0.37305748) q[0];
rz(2.864605) q[1];
sx q[1];
rz(-3.1412558) q[1];
sx q[1];
rz(-2.3854947) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33462229) q[0];
sx q[0];
rz(-0.96617132) q[0];
sx q[0];
rz(-2.0367677) q[0];
rz(-pi) q[1];
rz(0.71762849) q[2];
sx q[2];
rz(-0.52888479) q[2];
sx q[2];
rz(-2.8019049) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7623118) q[1];
sx q[1];
rz(-1.168615) q[1];
sx q[1];
rz(0.97907514) q[1];
x q[2];
rz(0.30310615) q[3];
sx q[3];
rz(-1.9744996) q[3];
sx q[3];
rz(-1.2991326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.21195352) q[2];
sx q[2];
rz(-0.55277199) q[2];
sx q[2];
rz(0.92145222) q[2];
rz(-0.067342162) q[3];
sx q[3];
rz(-1.9332956) q[3];
sx q[3];
rz(-1.4837846) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9901554) q[0];
sx q[0];
rz(-0.28011265) q[0];
sx q[0];
rz(-3.0113599) q[0];
rz(-0.79365927) q[1];
sx q[1];
rz(-3.1399813) q[1];
sx q[1];
rz(2.8614955) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5406427) q[0];
sx q[0];
rz(-1.4919992) q[0];
sx q[0];
rz(0.053683563) q[0];
rz(-0.21799223) q[2];
sx q[2];
rz(-1.6577066) q[2];
sx q[2];
rz(-0.72959585) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6895161) q[1];
sx q[1];
rz(-2.048564) q[1];
sx q[1];
rz(-2.6803451) q[1];
rz(-pi) q[2];
rz(-2.4305268) q[3];
sx q[3];
rz(-1.9424404) q[3];
sx q[3];
rz(1.4309331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.085999504) q[2];
sx q[2];
rz(-1.6722101) q[2];
sx q[2];
rz(-0.80297339) q[2];
rz(1.6064074) q[3];
sx q[3];
rz(-2.1740422) q[3];
sx q[3];
rz(-0.83139658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11237535) q[0];
sx q[0];
rz(-0.0028828415) q[0];
sx q[0];
rz(-0.10920864) q[0];
rz(0.373492) q[1];
sx q[1];
rz(-1.1871352) q[1];
sx q[1];
rz(-0.55140299) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5729061) q[0];
sx q[0];
rz(-0.99095067) q[0];
sx q[0];
rz(0.93007472) q[0];
x q[1];
rz(-0.18420561) q[2];
sx q[2];
rz(-1.6937243) q[2];
sx q[2];
rz(-0.21609989) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.046519) q[1];
sx q[1];
rz(-1.0380901) q[1];
sx q[1];
rz(-0.84360509) q[1];
x q[2];
rz(2.6104796) q[3];
sx q[3];
rz(-0.45980849) q[3];
sx q[3];
rz(-0.24548938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52770829) q[2];
sx q[2];
rz(-1.8324499) q[2];
sx q[2];
rz(-1.8180397) q[2];
rz(1.8771133) q[3];
sx q[3];
rz(-1.2885965) q[3];
sx q[3];
rz(0.0059787353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5453813) q[0];
sx q[0];
rz(-2.5061506) q[0];
sx q[0];
rz(2.4052461) q[0];
rz(0.21190724) q[1];
sx q[1];
rz(-2.1129463) q[1];
sx q[1];
rz(1.5440936) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6074267) q[0];
sx q[0];
rz(-1.5527927) q[0];
sx q[0];
rz(2.8624363) q[0];
x q[1];
rz(-0.030363516) q[2];
sx q[2];
rz(-1.5681475) q[2];
sx q[2];
rz(0.16636682) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5888583) q[1];
sx q[1];
rz(-0.97179669) q[1];
sx q[1];
rz(1.0042648) q[1];
rz(-pi) q[2];
rz(1.8964564) q[3];
sx q[3];
rz(-1.146637) q[3];
sx q[3];
rz(2.2617634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3046917) q[2];
sx q[2];
rz(-0.83536124) q[2];
sx q[2];
rz(-1.8677853) q[2];
rz(1.4443719) q[3];
sx q[3];
rz(-0.074051753) q[3];
sx q[3];
rz(1.6465638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.973751) q[0];
sx q[0];
rz(-1.583562) q[0];
sx q[0];
rz(-1.2927443) q[0];
rz(-1.5372859) q[1];
sx q[1];
rz(-2.2289386) q[1];
sx q[1];
rz(-2.9569721) q[1];
rz(0.040097728) q[2];
sx q[2];
rz(-1.5798777) q[2];
sx q[2];
rz(-0.29691534) q[2];
rz(0.96327412) q[3];
sx q[3];
rz(-1.6714851) q[3];
sx q[3];
rz(2.1906399) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
