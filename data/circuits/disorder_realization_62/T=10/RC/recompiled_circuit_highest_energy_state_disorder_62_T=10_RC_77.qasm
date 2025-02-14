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
rz(-0.57617968) q[0];
sx q[0];
rz(4.8186995) q[0];
sx q[0];
rz(8.7719593) q[0];
rz(2.6989812) q[1];
sx q[1];
rz(-1.1224597) q[1];
sx q[1];
rz(-2.519156) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1718423) q[0];
sx q[0];
rz(-1.6873708) q[0];
sx q[0];
rz(1.2571093) q[0];
rz(1.9046225) q[2];
sx q[2];
rz(-2.1458652) q[2];
sx q[2];
rz(1.0525557) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5099401) q[1];
sx q[1];
rz(-1.7544244) q[1];
sx q[1];
rz(-1.0650403) q[1];
rz(-pi) q[2];
rz(-0.962019) q[3];
sx q[3];
rz(-1.441027) q[3];
sx q[3];
rz(-2.6994841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6301253) q[2];
sx q[2];
rz(-1.3518535) q[2];
sx q[2];
rz(0.54568616) q[2];
rz(2.4557579) q[3];
sx q[3];
rz(-0.67982173) q[3];
sx q[3];
rz(-2.1849476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0551374) q[0];
sx q[0];
rz(-0.70239037) q[0];
sx q[0];
rz(1.0954683) q[0];
rz(-0.99758863) q[1];
sx q[1];
rz(-1.2350524) q[1];
sx q[1];
rz(1.197061) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2491978) q[0];
sx q[0];
rz(-2.4955242) q[0];
sx q[0];
rz(1.1419673) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95112339) q[2];
sx q[2];
rz(-2.24555) q[2];
sx q[2];
rz(-2.5017966) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6857288) q[1];
sx q[1];
rz(-1.2867754) q[1];
sx q[1];
rz(-1.0953254) q[1];
rz(-pi) q[2];
rz(1.7577111) q[3];
sx q[3];
rz(-1.2469562) q[3];
sx q[3];
rz(1.327654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7257488) q[2];
sx q[2];
rz(-1.7472605) q[2];
sx q[2];
rz(-1.3843298) q[2];
rz(-1.7978801) q[3];
sx q[3];
rz(-1.5735156) q[3];
sx q[3];
rz(1.2725007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62190732) q[0];
sx q[0];
rz(-2.3540731) q[0];
sx q[0];
rz(-2.7144879) q[0];
rz(-1.2318132) q[1];
sx q[1];
rz(-1.4991263) q[1];
sx q[1];
rz(-0.61418358) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1290951) q[0];
sx q[0];
rz(-1.3007727) q[0];
sx q[0];
rz(-2.8010023) q[0];
rz(0.74684192) q[2];
sx q[2];
rz(-1.421442) q[2];
sx q[2];
rz(0.17490444) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.28084785) q[1];
sx q[1];
rz(-0.13361803) q[1];
sx q[1];
rz(1.0722776) q[1];
x q[2];
rz(1.8133468) q[3];
sx q[3];
rz(-1.0772675) q[3];
sx q[3];
rz(2.7427615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.68982879) q[2];
sx q[2];
rz(-0.58388766) q[2];
sx q[2];
rz(0.27099398) q[2];
rz(-2.6594035) q[3];
sx q[3];
rz(-2.2995583) q[3];
sx q[3];
rz(1.2838001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.86376205) q[0];
sx q[0];
rz(-1.8067124) q[0];
sx q[0];
rz(2.7226287) q[0];
rz(2.9169182) q[1];
sx q[1];
rz(-1.2089665) q[1];
sx q[1];
rz(-3.0640501) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6200752) q[0];
sx q[0];
rz(-1.8635611) q[0];
sx q[0];
rz(-1.0241246) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9201905) q[2];
sx q[2];
rz(-2.3121693) q[2];
sx q[2];
rz(-1.0238613) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.78752199) q[1];
sx q[1];
rz(-1.7131299) q[1];
sx q[1];
rz(-2.700129) q[1];
rz(-pi) q[2];
rz(-0.61309149) q[3];
sx q[3];
rz(-0.35260751) q[3];
sx q[3];
rz(-2.4168454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.055723995) q[2];
sx q[2];
rz(-0.33129498) q[2];
sx q[2];
rz(1.9191939) q[2];
rz(-0.27082768) q[3];
sx q[3];
rz(-1.9041678) q[3];
sx q[3];
rz(2.6868611) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5643352) q[0];
sx q[0];
rz(-0.99522796) q[0];
sx q[0];
rz(-1.7876392) q[0];
rz(-2.1063781) q[1];
sx q[1];
rz(-1.926492) q[1];
sx q[1];
rz(1.6114906) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8299809) q[0];
sx q[0];
rz(-0.84048827) q[0];
sx q[0];
rz(0.82977022) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7351088) q[2];
sx q[2];
rz(-1.8949869) q[2];
sx q[2];
rz(2.6099082) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6275806) q[1];
sx q[1];
rz(-2.7360533) q[1];
sx q[1];
rz(-0.13579129) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6282998) q[3];
sx q[3];
rz(-1.1604476) q[3];
sx q[3];
rz(-1.8156605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3349541) q[2];
sx q[2];
rz(-2.1636476) q[2];
sx q[2];
rz(-0.063892603) q[2];
rz(-0.70819267) q[3];
sx q[3];
rz(-1.6876561) q[3];
sx q[3];
rz(-1.8074869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.719249) q[0];
sx q[0];
rz(-3.1184734) q[0];
sx q[0];
rz(-1.0759906) q[0];
rz(-3.0883582) q[1];
sx q[1];
rz(-0.85317555) q[1];
sx q[1];
rz(-2.7630973) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66010338) q[0];
sx q[0];
rz(-1.8253339) q[0];
sx q[0];
rz(1.4362364) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.041578023) q[2];
sx q[2];
rz(-1.6587306) q[2];
sx q[2];
rz(1.2961594) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.82847257) q[1];
sx q[1];
rz(-0.31658979) q[1];
sx q[1];
rz(0.60917882) q[1];
rz(2.1355576) q[3];
sx q[3];
rz(-0.88712403) q[3];
sx q[3];
rz(-2.5292252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1899015) q[2];
sx q[2];
rz(-1.4611763) q[2];
sx q[2];
rz(-2.107673) q[2];
rz(-2.2377491) q[3];
sx q[3];
rz(-2.240182) q[3];
sx q[3];
rz(-0.044053642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.0251625) q[0];
sx q[0];
rz(-1.6950386) q[0];
sx q[0];
rz(0.59301162) q[0];
rz(-3.016839) q[1];
sx q[1];
rz(-1.1589103) q[1];
sx q[1];
rz(-2.6483026) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.54087) q[0];
sx q[0];
rz(-2.713795) q[0];
sx q[0];
rz(-2.2912301) q[0];
rz(2.9134995) q[2];
sx q[2];
rz(-1.5950632) q[2];
sx q[2];
rz(0.14957854) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4932781) q[1];
sx q[1];
rz(-1.9946792) q[1];
sx q[1];
rz(-2.5878169) q[1];
rz(1.4629355) q[3];
sx q[3];
rz(-1.2965805) q[3];
sx q[3];
rz(1.8098243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.3844246) q[2];
sx q[2];
rz(-1.1606471) q[2];
sx q[2];
rz(1.9942795) q[2];
rz(0.82289639) q[3];
sx q[3];
rz(-1.1742679) q[3];
sx q[3];
rz(-0.66250044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5230474) q[0];
sx q[0];
rz(-1.86919) q[0];
sx q[0];
rz(-0.4185032) q[0];
rz(0.11323994) q[1];
sx q[1];
rz(-1.6721882) q[1];
sx q[1];
rz(1.07771) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1215873) q[0];
sx q[0];
rz(-1.5888056) q[0];
sx q[0];
rz(-1.2822582) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32025614) q[2];
sx q[2];
rz(-2.3204707) q[2];
sx q[2];
rz(2.1964354) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1220978) q[1];
sx q[1];
rz(-1.5217402) q[1];
sx q[1];
rz(8*pi/13) q[1];
rz(-pi) q[2];
rz(-1.5111132) q[3];
sx q[3];
rz(-2.2264997) q[3];
sx q[3];
rz(2.7969517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8354127) q[2];
sx q[2];
rz(-2.4330008) q[2];
sx q[2];
rz(-1.6861247) q[2];
rz(-0.16737394) q[3];
sx q[3];
rz(-0.51515976) q[3];
sx q[3];
rz(-1.0719871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48285943) q[0];
sx q[0];
rz(-1.8287683) q[0];
sx q[0];
rz(-0.40153781) q[0];
rz(-0.43831476) q[1];
sx q[1];
rz(-2.3500748) q[1];
sx q[1];
rz(-1.4567136) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0204937) q[0];
sx q[0];
rz(-1.5803845) q[0];
sx q[0];
rz(1.5367979) q[0];
rz(-1.3517604) q[2];
sx q[2];
rz(-0.92680762) q[2];
sx q[2];
rz(0.42962675) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5324133) q[1];
sx q[1];
rz(-1.0858156) q[1];
sx q[1];
rz(2.3663125) q[1];
x q[2];
rz(-2.0450122) q[3];
sx q[3];
rz(-2.0353531) q[3];
sx q[3];
rz(1.7024226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.38176408) q[2];
sx q[2];
rz(-0.20918736) q[2];
sx q[2];
rz(2.4299202) q[2];
rz(-3.107374) q[3];
sx q[3];
rz(-0.70845571) q[3];
sx q[3];
rz(-1.9862407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.3592247) q[0];
sx q[0];
rz(-1.490626) q[0];
sx q[0];
rz(2.3427298) q[0];
rz(-2.0955739) q[1];
sx q[1];
rz(-1.1963528) q[1];
sx q[1];
rz(0.23652133) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6470784) q[0];
sx q[0];
rz(-0.064726742) q[0];
sx q[0];
rz(-1.2184124) q[0];
rz(-pi) q[1];
rz(-2.4948214) q[2];
sx q[2];
rz(-2.1097217) q[2];
sx q[2];
rz(2.981271) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0933718) q[1];
sx q[1];
rz(-1.908675) q[1];
sx q[1];
rz(-0.45856234) q[1];
rz(3.0048851) q[3];
sx q[3];
rz(-2.4259704) q[3];
sx q[3];
rz(0.8807883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.11655) q[2];
sx q[2];
rz(-2.5849085) q[2];
sx q[2];
rz(0.64209783) q[2];
rz(-2.5721278) q[3];
sx q[3];
rz(-2.6537708) q[3];
sx q[3];
rz(0.08629442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4702598) q[0];
sx q[0];
rz(-1.3339806) q[0];
sx q[0];
rz(0.95680923) q[0];
rz(-1.9500465) q[1];
sx q[1];
rz(-1.5793431) q[1];
sx q[1];
rz(1.6021077) q[1];
rz(-1.446522) q[2];
sx q[2];
rz(-2.3588603) q[2];
sx q[2];
rz(-1.0609577) q[2];
rz(-1.996205) q[3];
sx q[3];
rz(-1.26401) q[3];
sx q[3];
rz(-1.1372529) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
