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
rz(0.25207818) q[0];
sx q[0];
rz(-2.5660388) q[0];
sx q[0];
rz(-0.12230305) q[0];
rz(-2.0343434) q[1];
sx q[1];
rz(-2.1802433) q[1];
sx q[1];
rz(0.038318757) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3275546) q[0];
sx q[0];
rz(-2.4310242) q[0];
sx q[0];
rz(2.1585805) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0228268) q[2];
sx q[2];
rz(-1.6890172) q[2];
sx q[2];
rz(-1.880065) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1520752) q[1];
sx q[1];
rz(-0.65836009) q[1];
sx q[1];
rz(1.0245628) q[1];
rz(3.0604912) q[3];
sx q[3];
rz(-0.24816324) q[3];
sx q[3];
rz(0.95141548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.64500874) q[2];
sx q[2];
rz(-1.4270447) q[2];
sx q[2];
rz(1.6928147) q[2];
rz(-0.19041666) q[3];
sx q[3];
rz(-1.0551635) q[3];
sx q[3];
rz(2.9922488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-1.9377624) q[0];
sx q[0];
rz(-1.8730524) q[0];
sx q[0];
rz(0.25800905) q[0];
rz(-1.5500655) q[1];
sx q[1];
rz(-1.9015046) q[1];
sx q[1];
rz(2.9749427) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7014871) q[0];
sx q[0];
rz(-0.99033725) q[0];
sx q[0];
rz(-1.6850182) q[0];
rz(-pi) q[1];
rz(0.3853674) q[2];
sx q[2];
rz(-1.3056985) q[2];
sx q[2];
rz(-2.6814658) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7609561) q[1];
sx q[1];
rz(-1.445862) q[1];
sx q[1];
rz(-2.4255468) q[1];
rz(-pi) q[2];
rz(1.0051425) q[3];
sx q[3];
rz(-2.0453302) q[3];
sx q[3];
rz(2.8090854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0698645) q[2];
sx q[2];
rz(-1.122415) q[2];
sx q[2];
rz(1.9094763) q[2];
rz(2.2169436) q[3];
sx q[3];
rz(-0.93622127) q[3];
sx q[3];
rz(-1.1054976) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.449618) q[0];
sx q[0];
rz(-1.469935) q[0];
sx q[0];
rz(-0.0019419226) q[0];
rz(-0.034189668) q[1];
sx q[1];
rz(-1.9296153) q[1];
sx q[1];
rz(-1.5984104) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1043248) q[0];
sx q[0];
rz(-1.1219949) q[0];
sx q[0];
rz(1.2842797) q[0];
rz(-0.43242792) q[2];
sx q[2];
rz(-0.78410599) q[2];
sx q[2];
rz(2.6570005) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.48422563) q[1];
sx q[1];
rz(-1.6099085) q[1];
sx q[1];
rz(2.7878321) q[1];
rz(-1.2528328) q[3];
sx q[3];
rz(-2.468716) q[3];
sx q[3];
rz(-1.8102136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89692846) q[2];
sx q[2];
rz(-2.7074773) q[2];
sx q[2];
rz(1.0373235) q[2];
rz(-2.0364929) q[3];
sx q[3];
rz(-1.5367855) q[3];
sx q[3];
rz(-1.1387811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26745519) q[0];
sx q[0];
rz(-0.48478165) q[0];
sx q[0];
rz(-1.7304035) q[0];
rz(-0.63938582) q[1];
sx q[1];
rz(-1.6849898) q[1];
sx q[1];
rz(3.0944518) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77098504) q[0];
sx q[0];
rz(-1.2895251) q[0];
sx q[0];
rz(-2.5384643) q[0];
rz(-pi) q[1];
rz(2.669692) q[2];
sx q[2];
rz(-1.7444897) q[2];
sx q[2];
rz(1.7431517) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85002181) q[1];
sx q[1];
rz(-0.98660368) q[1];
sx q[1];
rz(2.7136346) q[1];
rz(-pi) q[2];
rz(1.9023667) q[3];
sx q[3];
rz(-0.75115381) q[3];
sx q[3];
rz(2.0560665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3186657) q[2];
sx q[2];
rz(-1.1456127) q[2];
sx q[2];
rz(1.218943) q[2];
rz(-0.31560358) q[3];
sx q[3];
rz(-3.0867519) q[3];
sx q[3];
rz(-0.1489197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8197935) q[0];
sx q[0];
rz(-0.55235523) q[0];
sx q[0];
rz(-2.7929982) q[0];
rz(1.8888585) q[1];
sx q[1];
rz(-1.1628954) q[1];
sx q[1];
rz(-1.3016275) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4322223) q[0];
sx q[0];
rz(-0.8101058) q[0];
sx q[0];
rz(0.10929426) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5283777) q[2];
sx q[2];
rz(-2.4869031) q[2];
sx q[2];
rz(1.6729205) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.0685454) q[1];
sx q[1];
rz(-1.3497735) q[1];
sx q[1];
rz(-0.18494341) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0200649) q[3];
sx q[3];
rz(-1.2675084) q[3];
sx q[3];
rz(2.1316949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9584413) q[2];
sx q[2];
rz(-0.30854598) q[2];
sx q[2];
rz(-1.4754254) q[2];
rz(2.4557377) q[3];
sx q[3];
rz(-0.97969046) q[3];
sx q[3];
rz(-0.53387749) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1118065) q[0];
sx q[0];
rz(-0.15203467) q[0];
sx q[0];
rz(1.6492122) q[0];
rz(-0.97914186) q[1];
sx q[1];
rz(-1.6056332) q[1];
sx q[1];
rz(1.2243366) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5021073) q[0];
sx q[0];
rz(-1.7687135) q[0];
sx q[0];
rz(-1.5396126) q[0];
rz(-pi) q[1];
rz(1.3516828) q[2];
sx q[2];
rz(-2.1093528) q[2];
sx q[2];
rz(1.4736654) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.27367941) q[1];
sx q[1];
rz(-2.1845594) q[1];
sx q[1];
rz(3.1263208) q[1];
x q[2];
rz(1.1093418) q[3];
sx q[3];
rz(-0.68285817) q[3];
sx q[3];
rz(0.40753713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7225723) q[2];
sx q[2];
rz(-1.6442862) q[2];
sx q[2];
rz(-2.2255619) q[2];
rz(1.7570868) q[3];
sx q[3];
rz(-1.2576831) q[3];
sx q[3];
rz(0.70703435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0093805669) q[0];
sx q[0];
rz(-0.051932422) q[0];
sx q[0];
rz(-2.1873271) q[0];
rz(0.43680278) q[1];
sx q[1];
rz(-1.5635468) q[1];
sx q[1];
rz(-0.11016914) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5075275) q[0];
sx q[0];
rz(-1.332009) q[0];
sx q[0];
rz(2.8077543) q[0];
x q[1];
rz(-2.4442441) q[2];
sx q[2];
rz(-1.1432054) q[2];
sx q[2];
rz(-0.22874895) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4272144) q[1];
sx q[1];
rz(-2.7958779) q[1];
sx q[1];
rz(1.5835254) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0131936) q[3];
sx q[3];
rz(-0.54512776) q[3];
sx q[3];
rz(0.21827182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7097912) q[2];
sx q[2];
rz(-0.0063889901) q[2];
sx q[2];
rz(-0.94376454) q[2];
rz(-0.56378311) q[3];
sx q[3];
rz(-2.0457334) q[3];
sx q[3];
rz(-2.0311267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7875882) q[0];
sx q[0];
rz(-0.28689757) q[0];
sx q[0];
rz(0.34550825) q[0];
rz(1.0954789) q[1];
sx q[1];
rz(-1.4778719) q[1];
sx q[1];
rz(-1.5798205) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6248572) q[0];
sx q[0];
rz(-2.3586732) q[0];
sx q[0];
rz(1.5262414) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2279635) q[2];
sx q[2];
rz(-1.5217178) q[2];
sx q[2];
rz(0.14684453) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6553538) q[1];
sx q[1];
rz(-1.9353239) q[1];
sx q[1];
rz(-2.2097387) q[1];
x q[2];
rz(1.6027661) q[3];
sx q[3];
rz(-2.0684467) q[3];
sx q[3];
rz(-2.9835193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.76274189) q[2];
sx q[2];
rz(-2.8801535) q[2];
sx q[2];
rz(-1.4307107) q[2];
rz(0.53449574) q[3];
sx q[3];
rz(-1.4662687) q[3];
sx q[3];
rz(-0.9526332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6119824) q[0];
sx q[0];
rz(-0.069267608) q[0];
sx q[0];
rz(-1.8310504) q[0];
rz(2.2546841) q[1];
sx q[1];
rz(-2.3014018) q[1];
sx q[1];
rz(1.8720253) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5909605) q[0];
sx q[0];
rz(-2.7661588) q[0];
sx q[0];
rz(-1.9272789) q[0];
x q[1];
rz(-1.9475031) q[2];
sx q[2];
rz(-2.7393638) q[2];
sx q[2];
rz(0.62709432) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3514401) q[1];
sx q[1];
rz(-1.2987441) q[1];
sx q[1];
rz(2.9012009) q[1];
rz(-pi) q[2];
rz(-2.9379775) q[3];
sx q[3];
rz(-1.5316267) q[3];
sx q[3];
rz(2.7050582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5294007) q[2];
sx q[2];
rz(-1.8058913) q[2];
sx q[2];
rz(0.23294918) q[2];
rz(-1.8565146) q[3];
sx q[3];
rz(-2.464747) q[3];
sx q[3];
rz(-0.3860093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-1.9832298) q[0];
sx q[0];
rz(-0.60929275) q[0];
sx q[0];
rz(0.097271517) q[0];
rz(-0.80398792) q[1];
sx q[1];
rz(-0.709788) q[1];
sx q[1];
rz(-0.21673094) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6410206) q[0];
sx q[0];
rz(-1.4002698) q[0];
sx q[0];
rz(-3.0855623) q[0];
rz(-pi) q[1];
rz(2.8267536) q[2];
sx q[2];
rz(-0.14893571) q[2];
sx q[2];
rz(2.6268132) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.72783732) q[1];
sx q[1];
rz(-1.0168795) q[1];
sx q[1];
rz(-2.3604849) q[1];
x q[2];
rz(0.43301591) q[3];
sx q[3];
rz(-0.70954126) q[3];
sx q[3];
rz(-0.39329986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8486166) q[2];
sx q[2];
rz(-0.28102195) q[2];
sx q[2];
rz(-0.49523735) q[2];
rz(1.5589335) q[3];
sx q[3];
rz(-1.0844743) q[3];
sx q[3];
rz(-2.9787279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0611298) q[0];
sx q[0];
rz(-1.6292138) q[0];
sx q[0];
rz(-0.52393352) q[0];
rz(1.0864661) q[1];
sx q[1];
rz(-0.51849425) q[1];
sx q[1];
rz(0.35324221) q[1];
rz(-0.0072210014) q[2];
sx q[2];
rz(-1.5389969) q[2];
sx q[2];
rz(2.5015884) q[2];
rz(-0.85862715) q[3];
sx q[3];
rz(-1.5024019) q[3];
sx q[3];
rz(2.0988322) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
