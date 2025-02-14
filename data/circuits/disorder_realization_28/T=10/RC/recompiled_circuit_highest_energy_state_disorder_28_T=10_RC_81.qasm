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
rz(0.42143917) q[0];
sx q[0];
rz(5.0285664) q[0];
sx q[0];
rz(9.9528735) q[0];
rz(2.7948607) q[1];
sx q[1];
rz(-1.3277227) q[1];
sx q[1];
rz(0.2833856) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5534649) q[0];
sx q[0];
rz(-1.4736543) q[0];
sx q[0];
rz(-0.58421047) q[0];
rz(-0.82132116) q[2];
sx q[2];
rz(-1.9275713) q[2];
sx q[2];
rz(-2.5602464) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.99647616) q[1];
sx q[1];
rz(-1.3165425) q[1];
sx q[1];
rz(2.8338026) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3259726) q[3];
sx q[3];
rz(-2.7283165) q[3];
sx q[3];
rz(2.5101889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1957207) q[2];
sx q[2];
rz(-1.4737782) q[2];
sx q[2];
rz(0.26738581) q[2];
rz(-2.9684559) q[3];
sx q[3];
rz(-1.1370398) q[3];
sx q[3];
rz(0.6828298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.46083573) q[0];
sx q[0];
rz(-2.2956678) q[0];
sx q[0];
rz(-2.8435006) q[0];
rz(2.8744892) q[1];
sx q[1];
rz(-2.2593081) q[1];
sx q[1];
rz(2.2411236) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6783645) q[0];
sx q[0];
rz(-1.8186388) q[0];
sx q[0];
rz(2.8926587) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7212333) q[2];
sx q[2];
rz(-1.7152046) q[2];
sx q[2];
rz(0.23061801) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5067277) q[1];
sx q[1];
rz(-1.7657451) q[1];
sx q[1];
rz(2.4495803) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65761884) q[3];
sx q[3];
rz(-1.285835) q[3];
sx q[3];
rz(-2.5164925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4233826) q[2];
sx q[2];
rz(-0.98576468) q[2];
sx q[2];
rz(-0.59744376) q[2];
rz(-2.1256223) q[3];
sx q[3];
rz(-1.6705325) q[3];
sx q[3];
rz(1.5684675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51282561) q[0];
sx q[0];
rz(-0.91575423) q[0];
sx q[0];
rz(-3.0765007) q[0];
rz(2.0386631) q[1];
sx q[1];
rz(-1.256559) q[1];
sx q[1];
rz(-0.38806134) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93174879) q[0];
sx q[0];
rz(-1.8040621) q[0];
sx q[0];
rz(2.363028) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7746262) q[2];
sx q[2];
rz(-1.6410368) q[2];
sx q[2];
rz(-2.5052469) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1291099) q[1];
sx q[1];
rz(-1.7282244) q[1];
sx q[1];
rz(1.9858668) q[1];
x q[2];
rz(-2.2476419) q[3];
sx q[3];
rz(-2.7778447) q[3];
sx q[3];
rz(-1.8574024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1661561) q[2];
sx q[2];
rz(-1.0338975) q[2];
sx q[2];
rz(3.0530829) q[2];
rz(-2.4837808) q[3];
sx q[3];
rz(-2.7031873) q[3];
sx q[3];
rz(-1.2737761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45519644) q[0];
sx q[0];
rz(-0.96298591) q[0];
sx q[0];
rz(0.92940593) q[0];
rz(-1.9724253) q[1];
sx q[1];
rz(-0.87375748) q[1];
sx q[1];
rz(3.015231) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39462659) q[0];
sx q[0];
rz(-1.639059) q[0];
sx q[0];
rz(1.9771876) q[0];
rz(-2.8123145) q[2];
sx q[2];
rz(-1.3751404) q[2];
sx q[2];
rz(-0.31737993) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11495464) q[1];
sx q[1];
rz(-1.8787787) q[1];
sx q[1];
rz(2.6796209) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5764719) q[3];
sx q[3];
rz(-2.6905746) q[3];
sx q[3];
rz(-0.39263157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.29752842) q[2];
sx q[2];
rz(-1.0928096) q[2];
sx q[2];
rz(-1.7139942) q[2];
rz(1.1045688) q[3];
sx q[3];
rz(-1.5966871) q[3];
sx q[3];
rz(-0.67232084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70668689) q[0];
sx q[0];
rz(-0.66449419) q[0];
sx q[0];
rz(1.2427166) q[0];
rz(0.13936123) q[1];
sx q[1];
rz(-2.8262973) q[1];
sx q[1];
rz(2.349966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7032486) q[0];
sx q[0];
rz(-1.4171964) q[0];
sx q[0];
rz(3.1036882) q[0];
x q[1];
rz(1.7475732) q[2];
sx q[2];
rz(-1.0097754) q[2];
sx q[2];
rz(-1.068312) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0549889) q[1];
sx q[1];
rz(-1.5468925) q[1];
sx q[1];
rz(-1.8782645) q[1];
rz(-pi) q[2];
rz(-0.96782622) q[3];
sx q[3];
rz(-1.669291) q[3];
sx q[3];
rz(2.8525794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7255154) q[2];
sx q[2];
rz(-1.5308341) q[2];
sx q[2];
rz(-1.7804954) q[2];
rz(-1.037589) q[3];
sx q[3];
rz(-1.6890084) q[3];
sx q[3];
rz(3.1328372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1182275) q[0];
sx q[0];
rz(-2.8814377) q[0];
sx q[0];
rz(0.17876974) q[0];
rz(-2.692692) q[1];
sx q[1];
rz(-1.1652378) q[1];
sx q[1];
rz(1.9662205) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.878137) q[0];
sx q[0];
rz(-2.3317605) q[0];
sx q[0];
rz(1.4803588) q[0];
rz(-pi) q[1];
rz(2.3377445) q[2];
sx q[2];
rz(-1.7494252) q[2];
sx q[2];
rz(1.2527901) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7685942) q[1];
sx q[1];
rz(-1.914089) q[1];
sx q[1];
rz(2.3299135) q[1];
x q[2];
rz(0.39534335) q[3];
sx q[3];
rz(-2.5892777) q[3];
sx q[3];
rz(2.0022587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41345227) q[2];
sx q[2];
rz(-0.38273013) q[2];
sx q[2];
rz(0.21391301) q[2];
rz(1.4296069) q[3];
sx q[3];
rz(-1.471328) q[3];
sx q[3];
rz(-1.3235929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(1.0578617) q[0];
sx q[0];
rz(-1.2849176) q[0];
sx q[0];
rz(0.65823746) q[0];
rz(-1.4312875) q[1];
sx q[1];
rz(-0.88600102) q[1];
sx q[1];
rz(-1.5586982) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9607761) q[0];
sx q[0];
rz(-0.48191285) q[0];
sx q[0];
rz(-0.41367905) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8020242) q[2];
sx q[2];
rz(-1.2716846) q[2];
sx q[2];
rz(-2.0264152) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66872193) q[1];
sx q[1];
rz(-1.8897561) q[1];
sx q[1];
rz(-1.9265367) q[1];
x q[2];
rz(2.57929) q[3];
sx q[3];
rz(-0.65861693) q[3];
sx q[3];
rz(3.1059732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.3805286) q[2];
sx q[2];
rz(-1.5498127) q[2];
sx q[2];
rz(-1.2065678) q[2];
rz(-1.0208463) q[3];
sx q[3];
rz(-2.6267093) q[3];
sx q[3];
rz(2.2339039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9970488) q[0];
sx q[0];
rz(-2.9289065) q[0];
sx q[0];
rz(-0.3120684) q[0];
rz(2.8128305) q[1];
sx q[1];
rz(-1.6203251) q[1];
sx q[1];
rz(2.388248) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8341136) q[0];
sx q[0];
rz(-2.284482) q[0];
sx q[0];
rz(-0.77425787) q[0];
rz(-pi) q[1];
rz(-3.1198959) q[2];
sx q[2];
rz(-1.0602131) q[2];
sx q[2];
rz(-1.7006602) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.44239653) q[1];
sx q[1];
rz(-0.96328674) q[1];
sx q[1];
rz(2.1408666) q[1];
rz(2.3089319) q[3];
sx q[3];
rz(-0.7023905) q[3];
sx q[3];
rz(-2.663132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3030777) q[2];
sx q[2];
rz(-1.3884156) q[2];
sx q[2];
rz(-2.370749) q[2];
rz(-2.5924957) q[3];
sx q[3];
rz(-1.8531468) q[3];
sx q[3];
rz(0.87636605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86641208) q[0];
sx q[0];
rz(-2.488945) q[0];
sx q[0];
rz(-1.8021679) q[0];
rz(-2.0088947) q[1];
sx q[1];
rz(-0.39342543) q[1];
sx q[1];
rz(-0.045230953) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68494481) q[0];
sx q[0];
rz(-2.6620925) q[0];
sx q[0];
rz(1.9994237) q[0];
rz(-pi) q[1];
rz(-1.9275749) q[2];
sx q[2];
rz(-1.6851116) q[2];
sx q[2];
rz(-0.040093523) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9459234) q[1];
sx q[1];
rz(-1.821854) q[1];
sx q[1];
rz(1.4359479) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1858077) q[3];
sx q[3];
rz(-0.5866881) q[3];
sx q[3];
rz(-1.4101113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.905978) q[2];
sx q[2];
rz(-2.2048042) q[2];
sx q[2];
rz(-0.81602412) q[2];
rz(-1.7488545) q[3];
sx q[3];
rz(-2.0413028) q[3];
sx q[3];
rz(-1.5131153) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29869646) q[0];
sx q[0];
rz(-0.5492292) q[0];
sx q[0];
rz(1.5604875) q[0];
rz(2.9807978) q[1];
sx q[1];
rz(-2.000688) q[1];
sx q[1];
rz(-2.2198832) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7332382) q[0];
sx q[0];
rz(-1.7599765) q[0];
sx q[0];
rz(3.1319836) q[0];
rz(-2.3733487) q[2];
sx q[2];
rz(-2.3934954) q[2];
sx q[2];
rz(2.04271) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.28593496) q[1];
sx q[1];
rz(-1.287957) q[1];
sx q[1];
rz(0.15423473) q[1];
rz(-pi) q[2];
rz(0.62008786) q[3];
sx q[3];
rz(-0.26066565) q[3];
sx q[3];
rz(1.0083782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2738652) q[2];
sx q[2];
rz(-2.2160857) q[2];
sx q[2];
rz(1.072139) q[2];
rz(0.22465651) q[3];
sx q[3];
rz(-1.6499237) q[3];
sx q[3];
rz(2.2255161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0199725) q[0];
sx q[0];
rz(-2.2157123) q[0];
sx q[0];
rz(-2.8257688) q[0];
rz(-0.074180457) q[1];
sx q[1];
rz(-2.0769495) q[1];
sx q[1];
rz(3.1202797) q[1];
rz(-3.1138641) q[2];
sx q[2];
rz(-2.3578845) q[2];
sx q[2];
rz(3.1393928) q[2];
rz(2.737358) q[3];
sx q[3];
rz(-2.3334097) q[3];
sx q[3];
rz(2.4142053) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
