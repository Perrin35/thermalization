OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.83710837) q[0];
sx q[0];
rz(-1.4533071) q[0];
sx q[0];
rz(0.31153554) q[0];
rz(2.7040634) q[1];
sx q[1];
rz(-1.3181926) q[1];
sx q[1];
rz(-0.55895609) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5200978) q[0];
sx q[0];
rz(-1.1557475) q[0];
sx q[0];
rz(2.9893304) q[0];
rz(0.55732255) q[2];
sx q[2];
rz(-1.5601336) q[2];
sx q[2];
rz(2.9023841) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.47145876) q[1];
sx q[1];
rz(-2.11073) q[1];
sx q[1];
rz(-2.259841) q[1];
x q[2];
rz(-1.8564838) q[3];
sx q[3];
rz(-1.1477071) q[3];
sx q[3];
rz(0.74564122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3922334) q[2];
sx q[2];
rz(-1.2831251) q[2];
sx q[2];
rz(-0.63670811) q[2];
rz(2.2926245) q[3];
sx q[3];
rz(-2.520112) q[3];
sx q[3];
rz(0.23392114) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37671509) q[0];
sx q[0];
rz(-0.24704084) q[0];
sx q[0];
rz(-0.15287457) q[0];
rz(-2.3846467) q[1];
sx q[1];
rz(-1.5870973) q[1];
sx q[1];
rz(-0.98639948) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9372285) q[0];
sx q[0];
rz(-1.5110656) q[0];
sx q[0];
rz(1.6110957) q[0];
rz(-2.0979004) q[2];
sx q[2];
rz(-0.88000789) q[2];
sx q[2];
rz(1.8156798) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0806597) q[1];
sx q[1];
rz(-2.6903209) q[1];
sx q[1];
rz(0.56834759) q[1];
x q[2];
rz(0.13568474) q[3];
sx q[3];
rz(-1.2447345) q[3];
sx q[3];
rz(0.56860926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5622921) q[2];
sx q[2];
rz(-1.219517) q[2];
sx q[2];
rz(-2.3584649) q[2];
rz(0.018571818) q[3];
sx q[3];
rz(-1.637807) q[3];
sx q[3];
rz(2.7338681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.0884393) q[0];
sx q[0];
rz(-2.8494371) q[0];
sx q[0];
rz(-2.1799178) q[0];
rz(2.7812474) q[1];
sx q[1];
rz(-1.1018437) q[1];
sx q[1];
rz(0.12869421) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0592244) q[0];
sx q[0];
rz(-1.5229862) q[0];
sx q[0];
rz(1.6533018) q[0];
rz(2.1305069) q[2];
sx q[2];
rz(-2.2957544) q[2];
sx q[2];
rz(2.9618008) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3330831) q[1];
sx q[1];
rz(-1.2915478) q[1];
sx q[1];
rz(-0.46344325) q[1];
rz(-2.3388561) q[3];
sx q[3];
rz(-1.667913) q[3];
sx q[3];
rz(0.10955284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0814357) q[2];
sx q[2];
rz(-1.1971985) q[2];
sx q[2];
rz(-1.8998247) q[2];
rz(0.5870108) q[3];
sx q[3];
rz(-0.95654064) q[3];
sx q[3];
rz(-2.164042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63242763) q[0];
sx q[0];
rz(-2.2594663) q[0];
sx q[0];
rz(2.0571016) q[0];
rz(1.658461) q[1];
sx q[1];
rz(-0.56743923) q[1];
sx q[1];
rz(0.09253563) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4400892) q[0];
sx q[0];
rz(-1.2883696) q[0];
sx q[0];
rz(-2.861172) q[0];
rz(0.96524694) q[2];
sx q[2];
rz(-2.5778228) q[2];
sx q[2];
rz(-3.0639067) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1558518) q[1];
sx q[1];
rz(-2.7936613) q[1];
sx q[1];
rz(-2.9702529) q[1];
rz(-1.7792286) q[3];
sx q[3];
rz(-0.66550335) q[3];
sx q[3];
rz(0.044737577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3080421) q[2];
sx q[2];
rz(-1.7001067) q[2];
sx q[2];
rz(0.33205024) q[2];
rz(2.0856693) q[3];
sx q[3];
rz(-2.8639586) q[3];
sx q[3];
rz(2.5312996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32245359) q[0];
sx q[0];
rz(-1.2912913) q[0];
sx q[0];
rz(-0.21155587) q[0];
rz(1.3062723) q[1];
sx q[1];
rz(-1.897656) q[1];
sx q[1];
rz(0.64770118) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0487329) q[0];
sx q[0];
rz(-1.4591685) q[0];
sx q[0];
rz(-0.1603006) q[0];
rz(-pi) q[1];
rz(1.7251882) q[2];
sx q[2];
rz(-1.9152181) q[2];
sx q[2];
rz(1.7248578) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.16072893) q[1];
sx q[1];
rz(-0.54667066) q[1];
sx q[1];
rz(1.1854118) q[1];
rz(-pi) q[2];
rz(-1.5557489) q[3];
sx q[3];
rz(-1.5516557) q[3];
sx q[3];
rz(-1.0552989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.56090474) q[2];
sx q[2];
rz(-2.7320392) q[2];
sx q[2];
rz(-0.69331759) q[2];
rz(2.4723315) q[3];
sx q[3];
rz(-1.4368493) q[3];
sx q[3];
rz(-0.0049237331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82419056) q[0];
sx q[0];
rz(-2.9142002) q[0];
sx q[0];
rz(1.2325226) q[0];
rz(2.0690074) q[1];
sx q[1];
rz(-2.0697846) q[1];
sx q[1];
rz(0.17428621) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6264544) q[0];
sx q[0];
rz(-2.2720552) q[0];
sx q[0];
rz(-0.40909543) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3107489) q[2];
sx q[2];
rz(-2.2846662) q[2];
sx q[2];
rz(-1.6531528) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.78813997) q[1];
sx q[1];
rz(-0.8359682) q[1];
sx q[1];
rz(1.047903) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2586081) q[3];
sx q[3];
rz(-1.4600666) q[3];
sx q[3];
rz(-0.98715106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3198513) q[2];
sx q[2];
rz(-2.3345626) q[2];
sx q[2];
rz(-2.9439587) q[2];
rz(-2.8526784) q[3];
sx q[3];
rz(-0.87696004) q[3];
sx q[3];
rz(1.7355841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0777247) q[0];
sx q[0];
rz(-2.6517695) q[0];
sx q[0];
rz(-2.9329964) q[0];
rz(-2.1754307) q[1];
sx q[1];
rz(-1.9816793) q[1];
sx q[1];
rz(-1.5055515) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063232139) q[0];
sx q[0];
rz(-2.1713543) q[0];
sx q[0];
rz(-2.2230704) q[0];
rz(-3.0622919) q[2];
sx q[2];
rz(-2.7100025) q[2];
sx q[2];
rz(1.1232131) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5559616) q[1];
sx q[1];
rz(-2.0913843) q[1];
sx q[1];
rz(-0.52557892) q[1];
x q[2];
rz(-1.4671441) q[3];
sx q[3];
rz(-1.855587) q[3];
sx q[3];
rz(1.4485952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.85764) q[2];
sx q[2];
rz(-0.74308926) q[2];
sx q[2];
rz(-0.097578438) q[2];
rz(-1.7476667) q[3];
sx q[3];
rz(-1.7912309) q[3];
sx q[3];
rz(-0.71684366) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39032787) q[0];
sx q[0];
rz(-1.8126235) q[0];
sx q[0];
rz(0.51399291) q[0];
rz(-3.0184074) q[1];
sx q[1];
rz(-0.24736483) q[1];
sx q[1];
rz(0.93200144) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4829464) q[0];
sx q[0];
rz(-0.46447771) q[0];
sx q[0];
rz(-1.2645725) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2064662) q[2];
sx q[2];
rz(-1.7889651) q[2];
sx q[2];
rz(-0.93014923) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.276928) q[1];
sx q[1];
rz(-1.0229467) q[1];
sx q[1];
rz(2.1879556) q[1];
rz(0.19864948) q[3];
sx q[3];
rz(-2.2774787) q[3];
sx q[3];
rz(-2.228565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1121858) q[2];
sx q[2];
rz(-1.1108578) q[2];
sx q[2];
rz(-0.4294447) q[2];
rz(-1.9321692) q[3];
sx q[3];
rz(-2.7691787) q[3];
sx q[3];
rz(-2.4485574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76686239) q[0];
sx q[0];
rz(-1.6748036) q[0];
sx q[0];
rz(-1.8027579) q[0];
rz(-2.4354637) q[1];
sx q[1];
rz(-1.2495722) q[1];
sx q[1];
rz(-0.95058092) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75854077) q[0];
sx q[0];
rz(-2.0630815) q[0];
sx q[0];
rz(0.58347115) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.739903) q[2];
sx q[2];
rz(-0.42113129) q[2];
sx q[2];
rz(-0.087547628) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6420925) q[1];
sx q[1];
rz(-1.7576808) q[1];
sx q[1];
rz(1.0287813) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25456984) q[3];
sx q[3];
rz(-1.9707142) q[3];
sx q[3];
rz(-3.1018156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4799698) q[2];
sx q[2];
rz(-1.491549) q[2];
sx q[2];
rz(-2.6573112) q[2];
rz(-0.92710036) q[3];
sx q[3];
rz(-1.8323332) q[3];
sx q[3];
rz(-2.5203729) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9973688) q[0];
sx q[0];
rz(-0.08865083) q[0];
sx q[0];
rz(-0.22928672) q[0];
rz(-0.43481049) q[1];
sx q[1];
rz(-1.9186585) q[1];
sx q[1];
rz(2.4226709) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8573498) q[0];
sx q[0];
rz(-1.4765258) q[0];
sx q[0];
rz(-2.1988792) q[0];
rz(-pi) q[1];
rz(3.0232593) q[2];
sx q[2];
rz(-0.98630691) q[2];
sx q[2];
rz(0.10533939) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5303516) q[1];
sx q[1];
rz(-1.9396922) q[1];
sx q[1];
rz(1.7289274) q[1];
x q[2];
rz(-2.2273916) q[3];
sx q[3];
rz(-0.73540348) q[3];
sx q[3];
rz(-2.7587193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7583313) q[2];
sx q[2];
rz(-1.9394082) q[2];
sx q[2];
rz(-2.3948005) q[2];
rz(2.2693999) q[3];
sx q[3];
rz(-2.3132497) q[3];
sx q[3];
rz(-2.1993568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.223021) q[0];
sx q[0];
rz(-1.7445607) q[0];
sx q[0];
rz(1.8013409) q[0];
rz(-2.7643798) q[1];
sx q[1];
rz(-1.6709534) q[1];
sx q[1];
rz(0.51660641) q[1];
rz(0.5151667) q[2];
sx q[2];
rz(-2.5611521) q[2];
sx q[2];
rz(-0.44708154) q[2];
rz(-2.9711281) q[3];
sx q[3];
rz(-0.72453214) q[3];
sx q[3];
rz(1.7666045) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];