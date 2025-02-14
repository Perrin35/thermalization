OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2709687) q[0];
sx q[0];
rz(-0.55611098) q[0];
sx q[0];
rz(2.1882353) q[0];
rz(0.019729992) q[1];
sx q[1];
rz(-2.1058197) q[1];
sx q[1];
rz(-2.1929725) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23466388) q[0];
sx q[0];
rz(-1.0415823) q[0];
sx q[0];
rz(-1.0906903) q[0];
x q[1];
rz(-1.5679915) q[2];
sx q[2];
rz(-1.8082779) q[2];
sx q[2];
rz(-1.5495007) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.030101209) q[1];
sx q[1];
rz(-1.1426569) q[1];
sx q[1];
rz(-0.74049048) q[1];
rz(2.0998276) q[3];
sx q[3];
rz(-2.5210125) q[3];
sx q[3];
rz(2.4247501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3775776) q[2];
sx q[2];
rz(-0.074187584) q[2];
sx q[2];
rz(-2.3589676) q[2];
rz(-3.022656) q[3];
sx q[3];
rz(-2.1075893) q[3];
sx q[3];
rz(0.14757601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9453732) q[0];
sx q[0];
rz(-2.0202899) q[0];
sx q[0];
rz(-2.8837606) q[0];
rz(-0.091015426) q[1];
sx q[1];
rz(-1.0992522) q[1];
sx q[1];
rz(-1.6450504) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3548942) q[0];
sx q[0];
rz(-2.5845317) q[0];
sx q[0];
rz(-1.7003432) q[0];
rz(-pi) q[1];
rz(-0.74380959) q[2];
sx q[2];
rz(-0.48811661) q[2];
sx q[2];
rz(2.0632191) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1196668) q[1];
sx q[1];
rz(-2.7283784) q[1];
sx q[1];
rz(-1.6456804) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5037276) q[3];
sx q[3];
rz(-0.73299512) q[3];
sx q[3];
rz(0.021857787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0195134) q[2];
sx q[2];
rz(-1.8668819) q[2];
sx q[2];
rz(2.7369734) q[2];
rz(-0.98958611) q[3];
sx q[3];
rz(-1.8414958) q[3];
sx q[3];
rz(0.62304455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0290381) q[0];
sx q[0];
rz(-0.14415388) q[0];
sx q[0];
rz(2.3032904) q[0];
rz(-0.84838947) q[1];
sx q[1];
rz(-1.9219857) q[1];
sx q[1];
rz(-0.99260509) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.557195) q[0];
sx q[0];
rz(-0.8536754) q[0];
sx q[0];
rz(1.3846057) q[0];
rz(0.33432976) q[2];
sx q[2];
rz(-1.624776) q[2];
sx q[2];
rz(-2.8950429) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.015945399) q[1];
sx q[1];
rz(-1.9780668) q[1];
sx q[1];
rz(1.8512878) q[1];
x q[2];
rz(-0.29032536) q[3];
sx q[3];
rz(-1.6772406) q[3];
sx q[3];
rz(-0.31600472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8695716) q[2];
sx q[2];
rz(-1.973899) q[2];
sx q[2];
rz(-2.0167548) q[2];
rz(2.520842) q[3];
sx q[3];
rz(-2.1706332) q[3];
sx q[3];
rz(-0.11515215) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89897412) q[0];
sx q[0];
rz(-1.1599351) q[0];
sx q[0];
rz(0.17380357) q[0];
rz(0.68570343) q[1];
sx q[1];
rz(-1.655429) q[1];
sx q[1];
rz(-0.83820835) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80447557) q[0];
sx q[0];
rz(-0.67401471) q[0];
sx q[0];
rz(-2.3524257) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57212215) q[2];
sx q[2];
rz(-1.1231218) q[2];
sx q[2];
rz(-1.3037701) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.9058075) q[1];
sx q[1];
rz(-2.2407994) q[1];
sx q[1];
rz(1.7520105) q[1];
rz(-pi) q[2];
rz(-1.533314) q[3];
sx q[3];
rz(-2.0333383) q[3];
sx q[3];
rz(-1.9328062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.54245) q[2];
sx q[2];
rz(-1.4582381) q[2];
sx q[2];
rz(2.7899) q[2];
rz(1.9463978) q[3];
sx q[3];
rz(-1.9545133) q[3];
sx q[3];
rz(3.1174507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.46566063) q[0];
sx q[0];
rz(-2.6213578) q[0];
sx q[0];
rz(-2.7886673) q[0];
rz(0.30883166) q[1];
sx q[1];
rz(-2.1052723) q[1];
sx q[1];
rz(-3.1130863) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0519126) q[0];
sx q[0];
rz(-1.7453004) q[0];
sx q[0];
rz(-2.1580218) q[0];
x q[1];
rz(-0.79645313) q[2];
sx q[2];
rz(-2.6203794) q[2];
sx q[2];
rz(0.0040183684) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.18896477) q[1];
sx q[1];
rz(-0.82023925) q[1];
sx q[1];
rz(0.29839514) q[1];
rz(-2.0621544) q[3];
sx q[3];
rz(-0.48499987) q[3];
sx q[3];
rz(-0.62937832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9286524) q[2];
sx q[2];
rz(-3.0781367) q[2];
sx q[2];
rz(0.97079903) q[2];
rz(-1.0468696) q[3];
sx q[3];
rz(-1.4528843) q[3];
sx q[3];
rz(2.651732) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30188072) q[0];
sx q[0];
rz(-1.3588926) q[0];
sx q[0];
rz(3.0506328) q[0];
rz(1.92314) q[1];
sx q[1];
rz(-2.8107042) q[1];
sx q[1];
rz(-0.63922304) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6240182) q[0];
sx q[0];
rz(-1.7057452) q[0];
sx q[0];
rz(-0.3283511) q[0];
rz(-0.12219723) q[2];
sx q[2];
rz(-1.4336042) q[2];
sx q[2];
rz(2.4047763) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.51992937) q[1];
sx q[1];
rz(-0.93931336) q[1];
sx q[1];
rz(-2.2261621) q[1];
rz(-0.0019885824) q[3];
sx q[3];
rz(-1.244215) q[3];
sx q[3];
rz(2.0056564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0794534) q[2];
sx q[2];
rz(-2.1663351) q[2];
sx q[2];
rz(-2.6677168) q[2];
rz(-1.3300995) q[3];
sx q[3];
rz(-2.2606943) q[3];
sx q[3];
rz(2.7661095) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57629958) q[0];
sx q[0];
rz(-2.0602891) q[0];
sx q[0];
rz(0.22258776) q[0];
rz(1.0635771) q[1];
sx q[1];
rz(-1.5779326) q[1];
sx q[1];
rz(1.9482013) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.397307) q[0];
sx q[0];
rz(-0.82700506) q[0];
sx q[0];
rz(-1.4124407) q[0];
x q[1];
rz(-1.3435828) q[2];
sx q[2];
rz(-1.6889204) q[2];
sx q[2];
rz(2.3051777) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2031415) q[1];
sx q[1];
rz(-1.7060154) q[1];
sx q[1];
rz(-0.44288992) q[1];
rz(-pi) q[2];
rz(-0.16321793) q[3];
sx q[3];
rz(-1.3071408) q[3];
sx q[3];
rz(1.6632944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0003164) q[2];
sx q[2];
rz(-0.87339425) q[2];
sx q[2];
rz(2.5013962) q[2];
rz(0.021942465) q[3];
sx q[3];
rz(-1.7317737) q[3];
sx q[3];
rz(-0.35023165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6496277) q[0];
sx q[0];
rz(-1.3977298) q[0];
sx q[0];
rz(-2.6212027) q[0];
rz(0.87977663) q[1];
sx q[1];
rz(-1.2326515) q[1];
sx q[1];
rz(-0.89281503) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9111371) q[0];
sx q[0];
rz(-1.4833889) q[0];
sx q[0];
rz(1.5253517) q[0];
rz(-pi) q[1];
rz(-0.73142902) q[2];
sx q[2];
rz(-0.9584223) q[2];
sx q[2];
rz(-0.8548255) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.5984447) q[1];
sx q[1];
rz(-2.3259427) q[1];
sx q[1];
rz(-1.1786908) q[1];
x q[2];
rz(1.1475943) q[3];
sx q[3];
rz(-1.0931921) q[3];
sx q[3];
rz(1.3582317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2616547) q[2];
sx q[2];
rz(-2.0377906) q[2];
sx q[2];
rz(3.0547764) q[2];
rz(-0.9451198) q[3];
sx q[3];
rz(-2.7961531) q[3];
sx q[3];
rz(1.1300348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9574808) q[0];
sx q[0];
rz(-0.090395398) q[0];
sx q[0];
rz(-3.0775253) q[0];
rz(1.13569) q[1];
sx q[1];
rz(-1.4402729) q[1];
sx q[1];
rz(-0.51220977) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4762039) q[0];
sx q[0];
rz(-1.8955766) q[0];
sx q[0];
rz(0.43412125) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2733803) q[2];
sx q[2];
rz(-1.3033452) q[2];
sx q[2];
rz(-0.78154678) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.87742245) q[1];
sx q[1];
rz(-1.3119196) q[1];
sx q[1];
rz(-1.4822465) q[1];
rz(-1.5867611) q[3];
sx q[3];
rz(-1.9796238) q[3];
sx q[3];
rz(-0.4086993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3247437) q[2];
sx q[2];
rz(-2.4090359) q[2];
sx q[2];
rz(1.0653488) q[2];
rz(-1.6522853) q[3];
sx q[3];
rz(-0.95732147) q[3];
sx q[3];
rz(3.0491507) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5548993) q[0];
sx q[0];
rz(-0.41189343) q[0];
sx q[0];
rz(2.8107585) q[0];
rz(-1.2384442) q[1];
sx q[1];
rz(-0.60791433) q[1];
sx q[1];
rz(-0.27473658) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3392259) q[0];
sx q[0];
rz(-2.5792092) q[0];
sx q[0];
rz(2.2455162) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27336911) q[2];
sx q[2];
rz(-1.5393352) q[2];
sx q[2];
rz(-0.69845573) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3420136) q[1];
sx q[1];
rz(-1.430961) q[1];
sx q[1];
rz(-2.0796156) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4624743) q[3];
sx q[3];
rz(-0.59039718) q[3];
sx q[3];
rz(2.8670959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9773418) q[2];
sx q[2];
rz(-2.3873316) q[2];
sx q[2];
rz(1.3573307) q[2];
rz(2.6409798) q[3];
sx q[3];
rz(-1.5631792) q[3];
sx q[3];
rz(1.64465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(1.5951344) q[0];
sx q[0];
rz(-1.559579) q[0];
sx q[0];
rz(2.5478242) q[0];
rz(3.0733227) q[1];
sx q[1];
rz(-1.2208114) q[1];
sx q[1];
rz(0.023718871) q[1];
rz(-0.34577986) q[2];
sx q[2];
rz(-1.7263392) q[2];
sx q[2];
rz(1.6172258) q[2];
rz(-1.3832292) q[3];
sx q[3];
rz(-1.7456585) q[3];
sx q[3];
rz(-2.8432771) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
