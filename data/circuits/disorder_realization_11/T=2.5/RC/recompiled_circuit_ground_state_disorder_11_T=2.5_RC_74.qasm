OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1144855) q[0];
sx q[0];
rz(-0.067582421) q[0];
sx q[0];
rz(-0.58900589) q[0];
rz(-7.3569975) q[1];
sx q[1];
rz(4.8266657) q[1];
sx q[1];
rz(6.0756728) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02188259) q[0];
sx q[0];
rz(-1.8069022) q[0];
sx q[0];
rz(0.2050928) q[0];
x q[1];
rz(-1.4903388) q[2];
sx q[2];
rz(-1.3880037) q[2];
sx q[2];
rz(1.8535943) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7850951) q[1];
sx q[1];
rz(-1.5672088) q[1];
sx q[1];
rz(-1.5425578) q[1];
x q[2];
rz(-2.6785844) q[3];
sx q[3];
rz(-1.4742469) q[3];
sx q[3];
rz(0.90930401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.33240685) q[2];
sx q[2];
rz(-0.016533479) q[2];
sx q[2];
rz(-0.22745505) q[2];
rz(-2.4849232) q[3];
sx q[3];
rz(-0.69461099) q[3];
sx q[3];
rz(-0.53546661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6856132) q[0];
sx q[0];
rz(-0.032129012) q[0];
sx q[0];
rz(-2.6973714) q[0];
rz(2.0329068) q[1];
sx q[1];
rz(-1.8007092) q[1];
sx q[1];
rz(1.1681555) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21441169) q[0];
sx q[0];
rz(-2.4200038) q[0];
sx q[0];
rz(0.29436466) q[0];
rz(0.092165784) q[2];
sx q[2];
rz(-2.6861114) q[2];
sx q[2];
rz(-0.24126894) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.25623577) q[1];
sx q[1];
rz(-1.1215804) q[1];
sx q[1];
rz(-3.0577809) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2222667) q[3];
sx q[3];
rz(-0.72685234) q[3];
sx q[3];
rz(0.60692274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4593959) q[2];
sx q[2];
rz(-2.0708059) q[2];
sx q[2];
rz(-3.0219141) q[2];
rz(2.1648572) q[3];
sx q[3];
rz(-0.027438199) q[3];
sx q[3];
rz(-1.9500835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.479849) q[0];
sx q[0];
rz(-0.16525826) q[0];
sx q[0];
rz(1.4953493) q[0];
rz(-0.34548512) q[1];
sx q[1];
rz(-0.59436878) q[1];
sx q[1];
rz(-2.3631309) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.528196) q[0];
sx q[0];
rz(-1.5776618) q[0];
sx q[0];
rz(3.0860594) q[0];
x q[1];
rz(-0.81996452) q[2];
sx q[2];
rz(-1.1021309) q[2];
sx q[2];
rz(2.6033273) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.045588551) q[1];
sx q[1];
rz(-2.0109777) q[1];
sx q[1];
rz(0.11543302) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90980828) q[3];
sx q[3];
rz(-2.1902553) q[3];
sx q[3];
rz(-0.76515686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8567132) q[2];
sx q[2];
rz(-3.1046125) q[2];
sx q[2];
rz(-1.7516837) q[2];
rz(-0.022627929) q[3];
sx q[3];
rz(-0.32630625) q[3];
sx q[3];
rz(0.41445318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.49700272) q[0];
sx q[0];
rz(-3.1019326) q[0];
sx q[0];
rz(-0.52763754) q[0];
rz(0.35609326) q[1];
sx q[1];
rz(-1.5047319) q[1];
sx q[1];
rz(1.5632695) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2882345) q[0];
sx q[0];
rz(-0.70670595) q[0];
sx q[0];
rz(2.447829) q[0];
x q[1];
rz(-1.6354792) q[2];
sx q[2];
rz(-2.2475261) q[2];
sx q[2];
rz(-1.6874718) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9633052) q[1];
sx q[1];
rz(-1.7366341) q[1];
sx q[1];
rz(2.01009) q[1];
rz(3.0586518) q[3];
sx q[3];
rz(-2.0450052) q[3];
sx q[3];
rz(2.600973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.57033527) q[2];
sx q[2];
rz(-0.0084849914) q[2];
sx q[2];
rz(3.0961032) q[2];
rz(-0.82413834) q[3];
sx q[3];
rz(-2.6233311) q[3];
sx q[3];
rz(0.99388188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6028041) q[0];
sx q[0];
rz(-2.0863918) q[0];
sx q[0];
rz(-1.7428727) q[0];
rz(2.973373) q[1];
sx q[1];
rz(-1.8054447) q[1];
sx q[1];
rz(-0.22163637) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5580885) q[0];
sx q[0];
rz(-2.5595491) q[0];
sx q[0];
rz(-1.1314327) q[0];
rz(3.0732642) q[2];
sx q[2];
rz(-1.3573651) q[2];
sx q[2];
rz(2.1377856) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4532104) q[1];
sx q[1];
rz(-0.4391138) q[1];
sx q[1];
rz(-1.5136139) q[1];
x q[2];
rz(-0.92042376) q[3];
sx q[3];
rz(-1.2022965) q[3];
sx q[3];
rz(-0.96741048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.98442709) q[2];
sx q[2];
rz(-3.114279) q[2];
sx q[2];
rz(2.1402764) q[2];
rz(-2.2973513) q[3];
sx q[3];
rz(-3.0735569) q[3];
sx q[3];
rz(2.3943353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7897414) q[0];
sx q[0];
rz(-0.25775596) q[0];
sx q[0];
rz(-1.7716273) q[0];
rz(0.17579707) q[1];
sx q[1];
rz(-1.5134696) q[1];
sx q[1];
rz(-2.083185) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87230667) q[0];
sx q[0];
rz(-0.71505419) q[0];
sx q[0];
rz(1.1319517) q[0];
x q[1];
rz(2.1595512) q[2];
sx q[2];
rz(-1.5613404) q[2];
sx q[2];
rz(1.161425) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.26359162) q[1];
sx q[1];
rz(-1.8178838) q[1];
sx q[1];
rz(0.72551624) q[1];
x q[2];
rz(-2.2203683) q[3];
sx q[3];
rz(-1.8628253) q[3];
sx q[3];
rz(0.47367063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3093962) q[2];
sx q[2];
rz(-0.0038853566) q[2];
sx q[2];
rz(-0.84581172) q[2];
rz(-0.63515615) q[3];
sx q[3];
rz(-2.784415) q[3];
sx q[3];
rz(-0.14491189) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41088596) q[0];
sx q[0];
rz(-0.17179739) q[0];
sx q[0];
rz(-0.26096499) q[0];
rz(1.7102526) q[1];
sx q[1];
rz(-2.9899818) q[1];
sx q[1];
rz(1.3014911) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8943202) q[0];
sx q[0];
rz(-1.1218583) q[0];
sx q[0];
rz(-2.2821064) q[0];
rz(-pi) q[1];
rz(3.030483) q[2];
sx q[2];
rz(-1.5859563) q[2];
sx q[2];
rz(-0.76778256) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75682553) q[1];
sx q[1];
rz(-1.5315936) q[1];
sx q[1];
rz(0.43226166) q[1];
rz(0.00094125144) q[3];
sx q[3];
rz(-1.8347306) q[3];
sx q[3];
rz(-1.4838375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6053551) q[2];
sx q[2];
rz(-1.321512) q[2];
sx q[2];
rz(-0.028667299) q[2];
rz(3.0981433) q[3];
sx q[3];
rz(-2.9048007) q[3];
sx q[3];
rz(-1.6978327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0726149) q[0];
sx q[0];
rz(-0.34729877) q[0];
sx q[0];
rz(2.8436227) q[0];
rz(1.4120253) q[1];
sx q[1];
rz(-0.35930082) q[1];
sx q[1];
rz(1.5922458) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51247728) q[0];
sx q[0];
rz(-0.40226679) q[0];
sx q[0];
rz(-0.8921807) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.209433) q[2];
sx q[2];
rz(-3.123715) q[2];
sx q[2];
rz(1.0635384) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3181646) q[1];
sx q[1];
rz(-0.98814243) q[1];
sx q[1];
rz(1.2135439) q[1];
rz(-2.9202254) q[3];
sx q[3];
rz(-1.846617) q[3];
sx q[3];
rz(-2.3616977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90184244) q[2];
sx q[2];
rz(-3.1093137) q[2];
sx q[2];
rz(2.4816404) q[2];
rz(2.3615071) q[3];
sx q[3];
rz(-1.8476906) q[3];
sx q[3];
rz(-1.9735146) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50769794) q[0];
sx q[0];
rz(-0.036490353) q[0];
sx q[0];
rz(-0.80398917) q[0];
rz(-0.2313624) q[1];
sx q[1];
rz(-1.7295126) q[1];
sx q[1];
rz(3.0661327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7754525) q[0];
sx q[0];
rz(-2.7260927) q[0];
sx q[0];
rz(-0.52966161) q[0];
x q[1];
rz(-1.6310591) q[2];
sx q[2];
rz(-0.00060877006) q[2];
sx q[2];
rz(0.049545914) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4659404) q[1];
sx q[1];
rz(-2.2237101) q[1];
sx q[1];
rz(2.5616347) q[1];
rz(3.1072147) q[3];
sx q[3];
rz(-1.9581736) q[3];
sx q[3];
rz(-2.9971214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.017645322) q[2];
sx q[2];
rz(-0.46111527) q[2];
sx q[2];
rz(-2.7258415) q[2];
rz(1.6425284) q[3];
sx q[3];
rz(-0.00090229546) q[3];
sx q[3];
rz(2.3568962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7163664) q[0];
sx q[0];
rz(-0.056119053) q[0];
sx q[0];
rz(-0.29384336) q[0];
rz(1.5360688) q[1];
sx q[1];
rz(-2.2686281) q[1];
sx q[1];
rz(-1.6764838) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0241681) q[0];
sx q[0];
rz(-0.10483042) q[0];
sx q[0];
rz(-0.094314055) q[0];
x q[1];
rz(0.60428166) q[2];
sx q[2];
rz(-1.2245032) q[2];
sx q[2];
rz(-2.0228347) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0745213) q[1];
sx q[1];
rz(-1.5980729) q[1];
sx q[1];
rz(0.64665738) q[1];
rz(-1.5752389) q[3];
sx q[3];
rz(-1.5636139) q[3];
sx q[3];
rz(1.8069763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.49897033) q[2];
sx q[2];
rz(-0.62752807) q[2];
sx q[2];
rz(1.9210531) q[2];
rz(0.39310655) q[3];
sx q[3];
rz(-3.1246287) q[3];
sx q[3];
rz(0.19256798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6846234) q[0];
sx q[0];
rz(-1.4263117) q[0];
sx q[0];
rz(-0.34761467) q[0];
rz(-0.63812232) q[1];
sx q[1];
rz(-2.3633524) q[1];
sx q[1];
rz(-0.6755158) q[1];
rz(1.6498237) q[2];
sx q[2];
rz(-1.3776855) q[2];
sx q[2];
rz(2.9188271) q[2];
rz(-0.11278747) q[3];
sx q[3];
rz(-1.5776791) q[3];
sx q[3];
rz(0.37096544) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
