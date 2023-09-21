OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4797526) q[0];
sx q[0];
rz(-2.2979484) q[0];
sx q[0];
rz(2.9736829) q[0];
rz(1.1711988) q[1];
sx q[1];
rz(3.436915) q[1];
sx q[1];
rz(9.480939) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0060318) q[0];
sx q[0];
rz(-2.6132085) q[0];
sx q[0];
rz(-1.1392659) q[0];
rz(2.216823) q[2];
sx q[2];
rz(-1.7416818) q[2];
sx q[2];
rz(-2.8303763) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.72370428) q[1];
sx q[1];
rz(-1.0350409) q[1];
sx q[1];
rz(-2.8292311) q[1];
rz(-2.8817301) q[3];
sx q[3];
rz(-1.5279603) q[3];
sx q[3];
rz(-2.6916137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7636259) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(-2.7089233) q[2];
rz(-1.1928605) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(0.38309923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52779657) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(-1.8288076) q[0];
rz(-0.20547543) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(-1.1516494) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8173556) q[0];
sx q[0];
rz(-1.2756057) q[0];
sx q[0];
rz(-0.8582219) q[0];
rz(-pi) q[1];
rz(-2.0552911) q[2];
sx q[2];
rz(-2.0995579) q[2];
sx q[2];
rz(-2.8690086) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0484867) q[1];
sx q[1];
rz(-1.9890607) q[1];
sx q[1];
rz(-2.6622245) q[1];
x q[2];
rz(1.1851951) q[3];
sx q[3];
rz(-2.0413627) q[3];
sx q[3];
rz(-0.040599559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.0097222086) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(0.22182626) q[2];
rz(-2.7644073) q[3];
sx q[3];
rz(-0.42755869) q[3];
sx q[3];
rz(-2.3691573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31056988) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(-0.036852766) q[0];
rz(0.82551461) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(-0.056578606) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.409257) q[0];
sx q[0];
rz(-1.1261254) q[0];
sx q[0];
rz(2.1952941) q[0];
rz(-1.7043731) q[2];
sx q[2];
rz(-1.3165511) q[2];
sx q[2];
rz(0.1453407) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.11675662) q[1];
sx q[1];
rz(-2.7783238) q[1];
sx q[1];
rz(2.5519752) q[1];
rz(1.9583086) q[3];
sx q[3];
rz(-2.1790677) q[3];
sx q[3];
rz(1.1063948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8686707) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(-2.2154714) q[2];
rz(2.5849294) q[3];
sx q[3];
rz(-2.848048) q[3];
sx q[3];
rz(2.0986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91519231) q[0];
sx q[0];
rz(-0.72015786) q[0];
sx q[0];
rz(-2.3994989) q[0];
rz(1.1391976) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(-0.46359584) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6620561) q[0];
sx q[0];
rz(-2.2502796) q[0];
sx q[0];
rz(0.38328538) q[0];
rz(-pi) q[1];
rz(0.33275231) q[2];
sx q[2];
rz(-1.5126265) q[2];
sx q[2];
rz(1.0156877) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0193034) q[1];
sx q[1];
rz(-0.79320723) q[1];
sx q[1];
rz(1.3293468) q[1];
rz(-pi) q[2];
rz(-0.47834088) q[3];
sx q[3];
rz(-1.0676749) q[3];
sx q[3];
rz(-0.92418811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7745557) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(3.0920933) q[2];
rz(-0.1285304) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(-3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13609919) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(-0.29770011) q[0];
rz(0.4822576) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(-0.94435) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56086841) q[0];
sx q[0];
rz(-1.7797911) q[0];
sx q[0];
rz(0.046037721) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40446754) q[2];
sx q[2];
rz(-2.3239115) q[2];
sx q[2];
rz(0.58194619) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.22630616) q[1];
sx q[1];
rz(-1.5729135) q[1];
sx q[1];
rz(-1.6318984) q[1];
x q[2];
rz(1.827042) q[3];
sx q[3];
rz(-1.1503997) q[3];
sx q[3];
rz(0.30403578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.258761) q[2];
sx q[2];
rz(-1.0927039) q[2];
sx q[2];
rz(0.10822254) q[2];
rz(0.0023068874) q[3];
sx q[3];
rz(-1.6121515) q[3];
sx q[3];
rz(0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43679431) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(-2.5571402) q[0];
rz(-2.2553518) q[1];
sx q[1];
rz(-0.61683547) q[1];
sx q[1];
rz(3.086673) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14558218) q[0];
sx q[0];
rz(-1.8390363) q[0];
sx q[0];
rz(3.0560188) q[0];
rz(-pi) q[1];
rz(1.8940077) q[2];
sx q[2];
rz(-1.5530469) q[2];
sx q[2];
rz(-0.36349597) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0262895) q[1];
sx q[1];
rz(-2.4541828) q[1];
sx q[1];
rz(0.61647146) q[1];
x q[2];
rz(-0.92298569) q[3];
sx q[3];
rz(-0.69159782) q[3];
sx q[3];
rz(-1.5036316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3871258) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(1.0167271) q[2];
rz(-0.54404849) q[3];
sx q[3];
rz(-2.7816911) q[3];
sx q[3];
rz(1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.465437) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(-0.28453919) q[0];
rz(-2.1971205) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(-2.231266) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1658265) q[0];
sx q[0];
rz(-1.6253788) q[0];
sx q[0];
rz(-1.5058917) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1548642) q[2];
sx q[2];
rz(-1.8850733) q[2];
sx q[2];
rz(1.8679801) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6701723) q[1];
sx q[1];
rz(-0.52392611) q[1];
sx q[1];
rz(-0.36797221) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4207553) q[3];
sx q[3];
rz(-2.0257054) q[3];
sx q[3];
rz(-0.37973675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7515144) q[2];
sx q[2];
rz(-0.063515924) q[2];
sx q[2];
rz(0.92203036) q[2];
rz(-0.56728029) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(1.0197619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(-3.0122053) q[0];
rz(0.63240504) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(-0.30050373) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43113118) q[0];
sx q[0];
rz(-0.80695242) q[0];
sx q[0];
rz(1.0459082) q[0];
rz(0.76086107) q[2];
sx q[2];
rz(-2.0458474) q[2];
sx q[2];
rz(-2.8617815) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6943372) q[1];
sx q[1];
rz(-1.2303423) q[1];
sx q[1];
rz(1.7561595) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3281906) q[3];
sx q[3];
rz(-1.3374995) q[3];
sx q[3];
rz(-0.63026159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.58632103) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(-0.78197455) q[2];
rz(-2.590495) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(0.15792318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5732116) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(0.12776275) q[0];
rz(2.5993775) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(0.75884563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.471506) q[0];
sx q[0];
rz(-0.42450464) q[0];
sx q[0];
rz(1.5295117) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27300948) q[2];
sx q[2];
rz(-2.5148279) q[2];
sx q[2];
rz(-3.0221456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85941852) q[1];
sx q[1];
rz(-2.4255883) q[1];
sx q[1];
rz(-2.9031258) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5522478) q[3];
sx q[3];
rz(-2.4759001) q[3];
sx q[3];
rz(1.1902283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0163429) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(2.6515567) q[2];
rz(1.4222493) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816417) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(3.066257) q[0];
rz(-2.244859) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(2.5316701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41994914) q[0];
sx q[0];
rz(-2.3099265) q[0];
sx q[0];
rz(-1.182611) q[0];
rz(-pi) q[1];
rz(0.37656017) q[2];
sx q[2];
rz(-1.3137523) q[2];
sx q[2];
rz(3.0388447) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.48253179) q[1];
sx q[1];
rz(-0.7809124) q[1];
sx q[1];
rz(-0.34300967) q[1];
rz(-pi) q[2];
rz(-0.742357) q[3];
sx q[3];
rz(-1.1872113) q[3];
sx q[3];
rz(0.63391268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.909409) q[2];
sx q[2];
rz(-2.3164618) q[2];
sx q[2];
rz(-0.71371901) q[2];
rz(2.7632726) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(-0.87987125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80355766) q[0];
sx q[0];
rz(-1.9914347) q[0];
sx q[0];
rz(1.5557355) q[0];
rz(-0.65080416) q[1];
sx q[1];
rz(-1.4918068) q[1];
sx q[1];
rz(3.0202958) q[1];
rz(-0.030608721) q[2];
sx q[2];
rz(-1.3749214) q[2];
sx q[2];
rz(2.2236852) q[2];
rz(0.25272947) q[3];
sx q[3];
rz(-1.1128294) q[3];
sx q[3];
rz(2.2313234) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
