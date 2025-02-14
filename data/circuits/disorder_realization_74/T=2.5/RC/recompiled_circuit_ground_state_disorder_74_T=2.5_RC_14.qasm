OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.35535204) q[0];
sx q[0];
rz(-0.30682895) q[0];
sx q[0];
rz(-2.8428349) q[0];
rz(-0.6660676) q[1];
sx q[1];
rz(-0.63036418) q[1];
sx q[1];
rz(-1.6685553) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0399892) q[0];
sx q[0];
rz(-2.371006) q[0];
sx q[0];
rz(-1.6165074) q[0];
rz(-pi) q[1];
rz(-1.6285143) q[2];
sx q[2];
rz(-0.30122631) q[2];
sx q[2];
rz(1.1389552) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.3949497) q[1];
sx q[1];
rz(-2.267845) q[1];
sx q[1];
rz(-0.66956981) q[1];
rz(-pi) q[2];
rz(1.6410337) q[3];
sx q[3];
rz(-1.2956217) q[3];
sx q[3];
rz(-0.17632139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95639688) q[2];
sx q[2];
rz(-2.7608725) q[2];
sx q[2];
rz(-1.8308651) q[2];
rz(-3.0500566) q[3];
sx q[3];
rz(-2.5444578) q[3];
sx q[3];
rz(0.53102791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0700584) q[0];
sx q[0];
rz(-1.0907084) q[0];
sx q[0];
rz(-2.6614905) q[0];
rz(-1.9942888) q[1];
sx q[1];
rz(-0.6310178) q[1];
sx q[1];
rz(2.4400585) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.498256) q[0];
sx q[0];
rz(-1.0717046) q[0];
sx q[0];
rz(0.82478158) q[0];
rz(-pi) q[1];
rz(-0.84757324) q[2];
sx q[2];
rz(-0.59047359) q[2];
sx q[2];
rz(-2.5440885) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.50276668) q[1];
sx q[1];
rz(-2.354489) q[1];
sx q[1];
rz(0.9197286) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2626889) q[3];
sx q[3];
rz(-0.95617056) q[3];
sx q[3];
rz(-2.7648787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7445765) q[2];
sx q[2];
rz(-1.959356) q[2];
sx q[2];
rz(-1.5783295) q[2];
rz(0.27966106) q[3];
sx q[3];
rz(-2.4249488) q[3];
sx q[3];
rz(-2.7320812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79769832) q[0];
sx q[0];
rz(-0.41223031) q[0];
sx q[0];
rz(1.718234) q[0];
rz(-3.0176945) q[1];
sx q[1];
rz(-0.84404498) q[1];
sx q[1];
rz(-1.5843676) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14649571) q[0];
sx q[0];
rz(-2.2687758) q[0];
sx q[0];
rz(0.3905889) q[0];
rz(2.7742375) q[2];
sx q[2];
rz(-0.88161385) q[2];
sx q[2];
rz(0.85509461) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8516006) q[1];
sx q[1];
rz(-0.25892913) q[1];
sx q[1];
rz(-1.8845136) q[1];
rz(2.2267659) q[3];
sx q[3];
rz(-0.97712356) q[3];
sx q[3];
rz(-0.78730028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.065217) q[2];
sx q[2];
rz(-2.4880444) q[2];
sx q[2];
rz(-2.2067113) q[2];
rz(-2.8377418) q[3];
sx q[3];
rz(-0.6128208) q[3];
sx q[3];
rz(2.2936308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76871753) q[0];
sx q[0];
rz(-2.3621552) q[0];
sx q[0];
rz(-0.26707643) q[0];
rz(0.13485394) q[1];
sx q[1];
rz(-2.5649773) q[1];
sx q[1];
rz(-2.3042302) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98146472) q[0];
sx q[0];
rz(-0.61580393) q[0];
sx q[0];
rz(2.0439022) q[0];
rz(-pi) q[1];
rz(0.58394152) q[2];
sx q[2];
rz(-0.9204694) q[2];
sx q[2];
rz(-0.038829858) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3738651) q[1];
sx q[1];
rz(-0.30555913) q[1];
sx q[1];
rz(-0.69626804) q[1];
rz(-pi) q[2];
rz(1.6322953) q[3];
sx q[3];
rz(-1.5240747) q[3];
sx q[3];
rz(2.416147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5460633) q[2];
sx q[2];
rz(-0.17543051) q[2];
sx q[2];
rz(2.9626633) q[2];
rz(1.1692125) q[3];
sx q[3];
rz(-1.7763276) q[3];
sx q[3];
rz(-0.21807142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021407481) q[0];
sx q[0];
rz(-0.51161259) q[0];
sx q[0];
rz(2.2688493) q[0];
rz(1.1574289) q[1];
sx q[1];
rz(-0.87481421) q[1];
sx q[1];
rz(3.1246368) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30045045) q[0];
sx q[0];
rz(-0.73472856) q[0];
sx q[0];
rz(-0.018201306) q[0];
rz(-pi) q[1];
rz(-2.5418607) q[2];
sx q[2];
rz(-1.5625192) q[2];
sx q[2];
rz(-0.042009609) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.03093623) q[1];
sx q[1];
rz(-0.82582322) q[1];
sx q[1];
rz(0.40542545) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3153002) q[3];
sx q[3];
rz(-0.38681627) q[3];
sx q[3];
rz(1.7085962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.056328) q[2];
sx q[2];
rz(-1.7076098) q[2];
sx q[2];
rz(-2.8223574) q[2];
rz(0.6790092) q[3];
sx q[3];
rz(-1.346799) q[3];
sx q[3];
rz(2.3398248) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7769258) q[0];
sx q[0];
rz(-2.8209782) q[0];
sx q[0];
rz(-2.4989682) q[0];
rz(-1.369426) q[1];
sx q[1];
rz(-0.6232999) q[1];
sx q[1];
rz(-2.1646037) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8606023) q[0];
sx q[0];
rz(-0.15781584) q[0];
sx q[0];
rz(-0.67276038) q[0];
rz(-pi) q[1];
x q[1];
rz(2.634511) q[2];
sx q[2];
rz(-2.0424105) q[2];
sx q[2];
rz(-2.6952649) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7599267) q[1];
sx q[1];
rz(-1.5569512) q[1];
sx q[1];
rz(-0.87890505) q[1];
rz(-pi) q[2];
rz(-1.4897519) q[3];
sx q[3];
rz(-0.5872927) q[3];
sx q[3];
rz(2.0964266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6151108) q[2];
sx q[2];
rz(-1.7340163) q[2];
sx q[2];
rz(1.2283481) q[2];
rz(0.86722106) q[3];
sx q[3];
rz(-2.7286178) q[3];
sx q[3];
rz(-2.373608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7277471) q[0];
sx q[0];
rz(-2.9242046) q[0];
sx q[0];
rz(-0.17284285) q[0];
rz(-0.8051644) q[1];
sx q[1];
rz(-0.54904896) q[1];
sx q[1];
rz(0.13596143) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8589218) q[0];
sx q[0];
rz(-1.4513119) q[0];
sx q[0];
rz(-0.8651328) q[0];
rz(-pi) q[1];
rz(-1.682684) q[2];
sx q[2];
rz(-2.3896843) q[2];
sx q[2];
rz(-0.4900527) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31768979) q[1];
sx q[1];
rz(-3.025307) q[1];
sx q[1];
rz(-0.1648717) q[1];
rz(-pi) q[2];
x q[2];
rz(1.436266) q[3];
sx q[3];
rz(-1.8977652) q[3];
sx q[3];
rz(3.0643058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2148296) q[2];
sx q[2];
rz(-1.5462993) q[2];
sx q[2];
rz(0.1845486) q[2];
rz(-0.18937011) q[3];
sx q[3];
rz(-2.6034077) q[3];
sx q[3];
rz(-2.2034933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.1592584) q[0];
sx q[0];
rz(-0.33894798) q[0];
sx q[0];
rz(2.2151997) q[0];
rz(0.039208086) q[1];
sx q[1];
rz(-0.67752939) q[1];
sx q[1];
rz(2.1047986) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28473119) q[0];
sx q[0];
rz(-2.7960092) q[0];
sx q[0];
rz(2.1130108) q[0];
rz(-2.762297) q[2];
sx q[2];
rz(-0.46078983) q[2];
sx q[2];
rz(-2.5658105) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8421887) q[1];
sx q[1];
rz(-1.328598) q[1];
sx q[1];
rz(-0.41917015) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.751334) q[3];
sx q[3];
rz(-1.7446127) q[3];
sx q[3];
rz(2.9040608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.89256531) q[2];
sx q[2];
rz(-2.0165063) q[2];
sx q[2];
rz(0.79891515) q[2];
rz(-0.51698452) q[3];
sx q[3];
rz(-2.7345782) q[3];
sx q[3];
rz(-0.5504722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35140458) q[0];
sx q[0];
rz(-0.0049954448) q[0];
sx q[0];
rz(-1.5193526) q[0];
rz(1.5937357) q[1];
sx q[1];
rz(-0.82031119) q[1];
sx q[1];
rz(-0.69040745) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6192157) q[0];
sx q[0];
rz(-1.3478312) q[0];
sx q[0];
rz(-2.7821543) q[0];
rz(-pi) q[1];
x q[1];
rz(0.099748513) q[2];
sx q[2];
rz(-0.80173758) q[2];
sx q[2];
rz(-2.9005817) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.71296299) q[1];
sx q[1];
rz(-0.74364122) q[1];
sx q[1];
rz(2.3548467) q[1];
rz(-pi) q[2];
rz(1.810174) q[3];
sx q[3];
rz(-0.55174151) q[3];
sx q[3];
rz(2.5249979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4385628) q[2];
sx q[2];
rz(-2.5327693) q[2];
sx q[2];
rz(-2.1717066) q[2];
rz(0.25740933) q[3];
sx q[3];
rz(-2.9195547) q[3];
sx q[3];
rz(-1.7820057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53512204) q[0];
sx q[0];
rz(-2.4115925) q[0];
sx q[0];
rz(-0.085513376) q[0];
rz(-0.27587786) q[1];
sx q[1];
rz(-0.62092263) q[1];
sx q[1];
rz(1.1891018) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4258041) q[0];
sx q[0];
rz(-2.5809408) q[0];
sx q[0];
rz(0.7348357) q[0];
x q[1];
rz(1.5797473) q[2];
sx q[2];
rz(-1.8807966) q[2];
sx q[2];
rz(1.5218889) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.40285062) q[1];
sx q[1];
rz(-1.9704809) q[1];
sx q[1];
rz(0.12895361) q[1];
rz(0.34832921) q[3];
sx q[3];
rz(-1.4831721) q[3];
sx q[3];
rz(1.5548116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6659866) q[2];
sx q[2];
rz(-2.2624272) q[2];
sx q[2];
rz(-2.6574668) q[2];
rz(-0.13949805) q[3];
sx q[3];
rz(-2.3660584) q[3];
sx q[3];
rz(2.9805396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.465268) q[0];
sx q[0];
rz(-1.1558671) q[0];
sx q[0];
rz(-0.80094144) q[0];
rz(1.6839266) q[1];
sx q[1];
rz(-1.4977581) q[1];
sx q[1];
rz(-1.3118634) q[1];
rz(-0.64907907) q[2];
sx q[2];
rz(-1.4076283) q[2];
sx q[2];
rz(0.84244737) q[2];
rz(3.1112989) q[3];
sx q[3];
rz(-2.3843063) q[3];
sx q[3];
rz(1.6142308) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
