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
rz(-0.79688537) q[0];
sx q[0];
rz(-1.4112043) q[0];
sx q[0];
rz(-2.222173) q[0];
rz(1.3031651) q[1];
sx q[1];
rz(-0.04920955) q[1];
sx q[1];
rz(2.4278909) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.365276) q[0];
sx q[0];
rz(-2.8284024) q[0];
sx q[0];
rz(-2.0670783) q[0];
x q[1];
rz(-1.6993148) q[2];
sx q[2];
rz(-1.5084195) q[2];
sx q[2];
rz(-0.85278748) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9956869) q[1];
sx q[1];
rz(-2.6906081) q[1];
sx q[1];
rz(3.1383951) q[1];
rz(-0.034293745) q[3];
sx q[3];
rz(-2.0418797) q[3];
sx q[3];
rz(2.8054939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0539703) q[2];
sx q[2];
rz(-2.7555608) q[2];
sx q[2];
rz(-0.36960441) q[2];
rz(1.1450279) q[3];
sx q[3];
rz(-1.4729045) q[3];
sx q[3];
rz(-0.37231529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0358148) q[0];
sx q[0];
rz(-1.6749629) q[0];
sx q[0];
rz(-2.5685487) q[0];
rz(1.0967968) q[1];
sx q[1];
rz(-1.5410475) q[1];
sx q[1];
rz(0.47164741) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2315788) q[0];
sx q[0];
rz(-2.1699804) q[0];
sx q[0];
rz(-3.0082361) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4493982) q[2];
sx q[2];
rz(-1.0613777) q[2];
sx q[2];
rz(0.30864) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.009479) q[1];
sx q[1];
rz(-2.1188564) q[1];
sx q[1];
rz(2.4365042) q[1];
rz(-pi) q[2];
rz(-0.92949683) q[3];
sx q[3];
rz(-0.62493443) q[3];
sx q[3];
rz(2.9004824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7219217) q[2];
sx q[2];
rz(-2.2025043) q[2];
sx q[2];
rz(2.0527077) q[2];
rz(-1.0645083) q[3];
sx q[3];
rz(-1.1176611) q[3];
sx q[3];
rz(-2.815912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0629145) q[0];
sx q[0];
rz(-2.9525472) q[0];
sx q[0];
rz(-0.017070008) q[0];
rz(-0.71006376) q[1];
sx q[1];
rz(-0.83352572) q[1];
sx q[1];
rz(-0.56627083) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46895263) q[0];
sx q[0];
rz(-0.8297161) q[0];
sx q[0];
rz(-0.015182131) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.429661) q[2];
sx q[2];
rz(-1.1441137) q[2];
sx q[2];
rz(2.349145) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8108845) q[1];
sx q[1];
rz(-2.9363134) q[1];
sx q[1];
rz(2.7122871) q[1];
rz(-pi) q[2];
rz(0.26520437) q[3];
sx q[3];
rz(-1.6353161) q[3];
sx q[3];
rz(1.2539188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.337734) q[2];
sx q[2];
rz(-0.89581076) q[2];
sx q[2];
rz(0.0090948661) q[2];
rz(-3.005262) q[3];
sx q[3];
rz(-0.74458849) q[3];
sx q[3];
rz(1.2135308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.424054) q[0];
sx q[0];
rz(-0.67988765) q[0];
sx q[0];
rz(-2.9754382) q[0];
rz(-1.1135788) q[1];
sx q[1];
rz(-2.6515617) q[1];
sx q[1];
rz(0.16214935) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1095162) q[0];
sx q[0];
rz(-1.6639532) q[0];
sx q[0];
rz(-1.3901802) q[0];
x q[1];
rz(-1.5461363) q[2];
sx q[2];
rz(-1.0219821) q[2];
sx q[2];
rz(-1.8870712) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.6238193) q[1];
sx q[1];
rz(-2.2775067) q[1];
sx q[1];
rz(-1.0916187) q[1];
x q[2];
rz(-2.477263) q[3];
sx q[3];
rz(-1.2264612) q[3];
sx q[3];
rz(-0.20928247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98264155) q[2];
sx q[2];
rz(-1.0794159) q[2];
sx q[2];
rz(0.52784935) q[2];
rz(-2.2900901) q[3];
sx q[3];
rz(-2.8179171) q[3];
sx q[3];
rz(2.1984656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8714137) q[0];
sx q[0];
rz(-1.542792) q[0];
sx q[0];
rz(-2.340509) q[0];
rz(2.357645) q[1];
sx q[1];
rz(-0.74257094) q[1];
sx q[1];
rz(0.98091006) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0410515) q[0];
sx q[0];
rz(-1.9773736) q[0];
sx q[0];
rz(-1.340614) q[0];
x q[1];
rz(2.0774242) q[2];
sx q[2];
rz(-1.2384103) q[2];
sx q[2];
rz(1.4719964) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0438761) q[1];
sx q[1];
rz(-1.6642769) q[1];
sx q[1];
rz(3.082117) q[1];
rz(-2.122287) q[3];
sx q[3];
rz(-1.3634342) q[3];
sx q[3];
rz(2.1503748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0826147) q[2];
sx q[2];
rz(-1.5463983) q[2];
sx q[2];
rz(2.2583466) q[2];
rz(2.9412681) q[3];
sx q[3];
rz(-2.246558) q[3];
sx q[3];
rz(2.0506355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14466318) q[0];
sx q[0];
rz(-1.0522333) q[0];
sx q[0];
rz(-0.99217478) q[0];
rz(-0.80157533) q[1];
sx q[1];
rz(-1.8482607) q[1];
sx q[1];
rz(1.7209524) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29692263) q[0];
sx q[0];
rz(-1.2529182) q[0];
sx q[0];
rz(-0.65638377) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0367766) q[2];
sx q[2];
rz(-1.040852) q[2];
sx q[2];
rz(1.8304093) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.18773026) q[1];
sx q[1];
rz(-1.8988601) q[1];
sx q[1];
rz(0.66098722) q[1];
x q[2];
rz(1.6679974) q[3];
sx q[3];
rz(-1.5335113) q[3];
sx q[3];
rz(0.71747045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8617323) q[2];
sx q[2];
rz(-1.4272775) q[2];
sx q[2];
rz(0.53708616) q[2];
rz(3.1308657) q[3];
sx q[3];
rz(-0.74672943) q[3];
sx q[3];
rz(-0.85095325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2456197) q[0];
sx q[0];
rz(-2.8697822) q[0];
sx q[0];
rz(-0.33988345) q[0];
rz(1.8396359) q[1];
sx q[1];
rz(-2.1959031) q[1];
sx q[1];
rz(1.1713015) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.834496) q[0];
sx q[0];
rz(-2.8423873) q[0];
sx q[0];
rz(-0.12339573) q[0];
x q[1];
rz(0.76979678) q[2];
sx q[2];
rz(-2.8394921) q[2];
sx q[2];
rz(3.0506899) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75222941) q[1];
sx q[1];
rz(-3.1023355) q[1];
sx q[1];
rz(-2.1696349) q[1];
x q[2];
rz(0.24089916) q[3];
sx q[3];
rz(-2.5031075) q[3];
sx q[3];
rz(-1.0257667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12477144) q[2];
sx q[2];
rz(-1.7243959) q[2];
sx q[2];
rz(2.4701414) q[2];
rz(0.5365544) q[3];
sx q[3];
rz(-0.96009976) q[3];
sx q[3];
rz(-2.0898537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8272098) q[0];
sx q[0];
rz(-1.0541414) q[0];
sx q[0];
rz(-2.2453454) q[0];
rz(-0.25262901) q[1];
sx q[1];
rz(-2.2815506) q[1];
sx q[1];
rz(1.0505229) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4481215) q[0];
sx q[0];
rz(-0.52883178) q[0];
sx q[0];
rz(3.0419087) q[0];
rz(-pi) q[1];
rz(-2.4765313) q[2];
sx q[2];
rz(-1.6740693) q[2];
sx q[2];
rz(-2.2879083) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3747762) q[1];
sx q[1];
rz(-1.720806) q[1];
sx q[1];
rz(2.0376202) q[1];
rz(2.3837621) q[3];
sx q[3];
rz(-1.5327454) q[3];
sx q[3];
rz(-1.1421957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5184021) q[2];
sx q[2];
rz(-2.9727327) q[2];
sx q[2];
rz(1.5625578) q[2];
rz(1.3671499) q[3];
sx q[3];
rz(-2.0953777) q[3];
sx q[3];
rz(-2.3794543) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6147181) q[0];
sx q[0];
rz(-1.7521097) q[0];
sx q[0];
rz(-2.7516464) q[0];
rz(-0.82603106) q[1];
sx q[1];
rz(-2.4470058) q[1];
sx q[1];
rz(1.0521851) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56985939) q[0];
sx q[0];
rz(-2.2433407) q[0];
sx q[0];
rz(2.74617) q[0];
x q[1];
rz(2.7446943) q[2];
sx q[2];
rz(-1.7328615) q[2];
sx q[2];
rz(1.6610749) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.70665592) q[1];
sx q[1];
rz(-1.860992) q[1];
sx q[1];
rz(-1.9761843) q[1];
x q[2];
rz(0.47635079) q[3];
sx q[3];
rz(-2.1423376) q[3];
sx q[3];
rz(-1.0494572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0612001) q[2];
sx q[2];
rz(-1.281851) q[2];
sx q[2];
rz(-2.6123987) q[2];
rz(2.3906294) q[3];
sx q[3];
rz(-0.96143985) q[3];
sx q[3];
rz(-2.9474337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.879409) q[0];
sx q[0];
rz(-0.59469596) q[0];
sx q[0];
rz(1.9770812) q[0];
rz(1.2871845) q[1];
sx q[1];
rz(-2.4559805) q[1];
sx q[1];
rz(-0.50416344) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7499649) q[0];
sx q[0];
rz(-0.34191439) q[0];
sx q[0];
rz(2.2750399) q[0];
rz(-1.590056) q[2];
sx q[2];
rz(-1.9474721) q[2];
sx q[2];
rz(-0.45805675) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.530624) q[1];
sx q[1];
rz(-1.8357616) q[1];
sx q[1];
rz(0.19719657) q[1];
rz(-0.88621288) q[3];
sx q[3];
rz(-1.7401812) q[3];
sx q[3];
rz(-2.2365295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9265499) q[2];
sx q[2];
rz(-1.5861009) q[2];
sx q[2];
rz(-0.034493383) q[2];
rz(1.5283594) q[3];
sx q[3];
rz(-2.4210763) q[3];
sx q[3];
rz(-0.31149402) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16019776) q[0];
sx q[0];
rz(-1.5481411) q[0];
sx q[0];
rz(2.7465469) q[0];
rz(2.1389217) q[1];
sx q[1];
rz(-1.6802588) q[1];
sx q[1];
rz(-0.37793876) q[1];
rz(2.6832081) q[2];
sx q[2];
rz(-1.7469363) q[2];
sx q[2];
rz(1.7652702) q[2];
rz(0.18319753) q[3];
sx q[3];
rz(-0.33231322) q[3];
sx q[3];
rz(-0.65381369) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
