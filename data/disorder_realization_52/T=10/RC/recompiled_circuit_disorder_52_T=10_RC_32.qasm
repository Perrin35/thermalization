OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5157226) q[0];
sx q[0];
rz(-0.54870257) q[0];
sx q[0];
rz(-0.8843511) q[0];
rz(-1.7110775) q[1];
sx q[1];
rz(-0.95354748) q[1];
sx q[1];
rz(1.6391099) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26412548) q[0];
sx q[0];
rz(-1.2116417) q[0];
sx q[0];
rz(-0.19192625) q[0];
rz(-3.1332364) q[2];
sx q[2];
rz(-0.5651606) q[2];
sx q[2];
rz(-1.9985808) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1391746) q[1];
sx q[1];
rz(-2.8867509) q[1];
sx q[1];
rz(0.26267085) q[1];
rz(-pi) q[2];
rz(0.35688551) q[3];
sx q[3];
rz(-2.2483453) q[3];
sx q[3];
rz(2.4524636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.78645906) q[2];
sx q[2];
rz(-2.3278475) q[2];
sx q[2];
rz(2.4856429) q[2];
rz(1.2077228) q[3];
sx q[3];
rz(-1.9658807) q[3];
sx q[3];
rz(2.1470127) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2475964) q[0];
sx q[0];
rz(-2.7212454) q[0];
sx q[0];
rz(-0.43352747) q[0];
rz(2.9128089) q[1];
sx q[1];
rz(-0.42963916) q[1];
sx q[1];
rz(3.1343592) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9418966) q[0];
sx q[0];
rz(-1.4775839) q[0];
sx q[0];
rz(1.9641563) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3184001) q[2];
sx q[2];
rz(-2.7232812) q[2];
sx q[2];
rz(2.7559466) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0011598) q[1];
sx q[1];
rz(-1.7908485) q[1];
sx q[1];
rz(0.37186719) q[1];
rz(-0.98874436) q[3];
sx q[3];
rz(-2.3855004) q[3];
sx q[3];
rz(-2.3088629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5269512) q[2];
sx q[2];
rz(-0.80792892) q[2];
sx q[2];
rz(2.4439404) q[2];
rz(0.12156045) q[3];
sx q[3];
rz(-1.2391042) q[3];
sx q[3];
rz(0.30383032) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4678629) q[0];
sx q[0];
rz(-2.2255852) q[0];
sx q[0];
rz(1.3695705) q[0];
rz(1.2415775) q[1];
sx q[1];
rz(-1.7280271) q[1];
sx q[1];
rz(0.26161584) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81944377) q[0];
sx q[0];
rz(-1.5804287) q[0];
sx q[0];
rz(-1.5887512) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82654731) q[2];
sx q[2];
rz(-1.6625704) q[2];
sx q[2];
rz(-1.9020136) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97754543) q[1];
sx q[1];
rz(-2.3512212) q[1];
sx q[1];
rz(-0.17069823) q[1];
x q[2];
rz(-1.6022127) q[3];
sx q[3];
rz(-1.0580225) q[3];
sx q[3];
rz(-2.153742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4613142) q[2];
sx q[2];
rz(-1.6363982) q[2];
sx q[2];
rz(0.20351163) q[2];
rz(0.92173785) q[3];
sx q[3];
rz(-1.8739871) q[3];
sx q[3];
rz(-0.27954277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1638284) q[0];
sx q[0];
rz(-1.5638567) q[0];
sx q[0];
rz(-1.5699566) q[0];
rz(2.1381901) q[1];
sx q[1];
rz(-1.827821) q[1];
sx q[1];
rz(-1.2483695) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0700064) q[0];
sx q[0];
rz(-1.4741815) q[0];
sx q[0];
rz(-3.0759401) q[0];
rz(-1.3893045) q[2];
sx q[2];
rz(-1.3647807) q[2];
sx q[2];
rz(-0.69603053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.26989386) q[1];
sx q[1];
rz(-1.6819685) q[1];
sx q[1];
rz(-0.65472366) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8321886) q[3];
sx q[3];
rz(-2.2769615) q[3];
sx q[3];
rz(-0.9725001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0288329) q[2];
sx q[2];
rz(-1.8303822) q[2];
sx q[2];
rz(-2.3045585) q[2];
rz(1.933243) q[3];
sx q[3];
rz(-1.2669867) q[3];
sx q[3];
rz(2.3560431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1355302) q[0];
sx q[0];
rz(-1.0536138) q[0];
sx q[0];
rz(-0.77520448) q[0];
rz(-0.40183055) q[1];
sx q[1];
rz(-2.1907175) q[1];
sx q[1];
rz(0.90243375) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1989312) q[0];
sx q[0];
rz(-0.64093243) q[0];
sx q[0];
rz(3.065227) q[0];
rz(-pi) q[1];
rz(-0.4544223) q[2];
sx q[2];
rz(-1.7917969) q[2];
sx q[2];
rz(-2.1511252) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.93859766) q[1];
sx q[1];
rz(-1.0293759) q[1];
sx q[1];
rz(-0.90403892) q[1];
rz(-1.5794472) q[3];
sx q[3];
rz(-2.8248565) q[3];
sx q[3];
rz(-0.99416152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9465785) q[2];
sx q[2];
rz(-2.6533551) q[2];
sx q[2];
rz(1.9449332) q[2];
rz(-1.6992016) q[3];
sx q[3];
rz(-1.8701575) q[3];
sx q[3];
rz(1.6825914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8114132) q[0];
sx q[0];
rz(-2.9646962) q[0];
sx q[0];
rz(-0.50317558) q[0];
rz(1.6852089) q[1];
sx q[1];
rz(-2.0676985) q[1];
sx q[1];
rz(2.9398289) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0148894) q[0];
sx q[0];
rz(-3.0198583) q[0];
sx q[0];
rz(0.99862167) q[0];
x q[1];
rz(2.6642338) q[2];
sx q[2];
rz(-1.1672033) q[2];
sx q[2];
rz(-1.8237643) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6946053) q[1];
sx q[1];
rz(-1.1519377) q[1];
sx q[1];
rz(2.6161731) q[1];
rz(-pi) q[2];
rz(-1.1062578) q[3];
sx q[3];
rz(-0.89336508) q[3];
sx q[3];
rz(-0.57074742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9266944) q[2];
sx q[2];
rz(-1.5025257) q[2];
sx q[2];
rz(-2.8743437) q[2];
rz(0.8231419) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(0.87583035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8478407) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(-3.0840432) q[0];
rz(1.6607704) q[1];
sx q[1];
rz(-1.7819504) q[1];
sx q[1];
rz(2.1988791) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38521117) q[0];
sx q[0];
rz(-2.4009973) q[0];
sx q[0];
rz(-2.5368607) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16072388) q[2];
sx q[2];
rz(-1.7956927) q[2];
sx q[2];
rz(-2.7149534) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.47989935) q[1];
sx q[1];
rz(-2.6091895) q[1];
sx q[1];
rz(-0.9049306) q[1];
rz(-1.1293291) q[3];
sx q[3];
rz(-0.55674508) q[3];
sx q[3];
rz(2.7260775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0223579) q[2];
sx q[2];
rz(-1.0417577) q[2];
sx q[2];
rz(-2.7589202) q[2];
rz(-1.0391957) q[3];
sx q[3];
rz(-1.8363876) q[3];
sx q[3];
rz(2.1634845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7169749) q[0];
sx q[0];
rz(-3.0452403) q[0];
sx q[0];
rz(-2.8714645) q[0];
rz(2.5121636) q[1];
sx q[1];
rz(-2.4286178) q[1];
sx q[1];
rz(0.28392917) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1536381) q[0];
sx q[0];
rz(-2.0001912) q[0];
sx q[0];
rz(0.9149787) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21247272) q[2];
sx q[2];
rz(-0.61056911) q[2];
sx q[2];
rz(-3.0688822) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.64360196) q[1];
sx q[1];
rz(-2.0703348) q[1];
sx q[1];
rz(1.040578) q[1];
x q[2];
rz(-0.038589434) q[3];
sx q[3];
rz(-0.66017294) q[3];
sx q[3];
rz(0.71036464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55390629) q[2];
sx q[2];
rz(-0.17051414) q[2];
sx q[2];
rz(1.930687) q[2];
rz(-2.8816913) q[3];
sx q[3];
rz(-0.62429684) q[3];
sx q[3];
rz(-0.62817812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0652086) q[0];
sx q[0];
rz(-2.5807091) q[0];
sx q[0];
rz(-0.22924766) q[0];
rz(-0.30300888) q[1];
sx q[1];
rz(-1.3906994) q[1];
sx q[1];
rz(-1.4607666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54366771) q[0];
sx q[0];
rz(-0.28408465) q[0];
sx q[0];
rz(3.0790867) q[0];
rz(0.37947189) q[2];
sx q[2];
rz(-1.737829) q[2];
sx q[2];
rz(-0.52969474) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9257099) q[1];
sx q[1];
rz(-2.6162418) q[1];
sx q[1];
rz(-3.0970296) q[1];
rz(-pi) q[2];
rz(-1.1838412) q[3];
sx q[3];
rz(-2.317252) q[3];
sx q[3];
rz(1.6798306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4298657) q[2];
sx q[2];
rz(-1.9248328) q[2];
sx q[2];
rz(2.4460068) q[2];
rz(-2.7097278) q[3];
sx q[3];
rz(-0.46468195) q[3];
sx q[3];
rz(0.7152043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5678976) q[0];
sx q[0];
rz(-1.9298113) q[0];
sx q[0];
rz(2.8046872) q[0];
rz(2.9341872) q[1];
sx q[1];
rz(-2.1131056) q[1];
sx q[1];
rz(-2.7609603) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99285179) q[0];
sx q[0];
rz(-1.7721575) q[0];
sx q[0];
rz(-3.0701748) q[0];
x q[1];
rz(0.74110465) q[2];
sx q[2];
rz(-1.1198824) q[2];
sx q[2];
rz(-1.8795183) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29585782) q[1];
sx q[1];
rz(-1.0746135) q[1];
sx q[1];
rz(-0.22175281) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5467039) q[3];
sx q[3];
rz(-2.4747362) q[3];
sx q[3];
rz(2.8951333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1404861) q[2];
sx q[2];
rz(-1.9753186) q[2];
sx q[2];
rz(-2.005119) q[2];
rz(-3.100637) q[3];
sx q[3];
rz(-0.80871964) q[3];
sx q[3];
rz(-1.827318) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8052335) q[0];
sx q[0];
rz(-2.2362066) q[0];
sx q[0];
rz(2.8295828) q[0];
rz(-2.1144755) q[1];
sx q[1];
rz(-1.849091) q[1];
sx q[1];
rz(-1.0277933) q[1];
rz(2.8935511) q[2];
sx q[2];
rz(-1.7938062) q[2];
sx q[2];
rz(1.2563406) q[2];
rz(1.9789226) q[3];
sx q[3];
rz(-0.47331953) q[3];
sx q[3];
rz(2.4389653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
