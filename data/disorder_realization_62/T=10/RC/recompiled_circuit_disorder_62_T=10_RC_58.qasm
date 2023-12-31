OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9010889) q[0];
sx q[0];
rz(-0.17740372) q[0];
sx q[0];
rz(1.1343962) q[0];
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(0.66361767) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21270574) q[0];
sx q[0];
rz(-1.1097254) q[0];
sx q[0];
rz(-1.5972208) q[0];
rz(-0.77779777) q[2];
sx q[2];
rz(-1.8282837) q[2];
sx q[2];
rz(-2.4553026) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9259778) q[1];
sx q[1];
rz(-1.3091334) q[1];
sx q[1];
rz(1.2222626) q[1];
x q[2];
rz(1.1575559) q[3];
sx q[3];
rz(-2.8765656) q[3];
sx q[3];
rz(-2.7812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.87876451) q[2];
sx q[2];
rz(-2.6972013) q[2];
sx q[2];
rz(-3.0905241) q[2];
rz(2.5845394) q[3];
sx q[3];
rz(-0.80018187) q[3];
sx q[3];
rz(1.5867656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5550845) q[0];
sx q[0];
rz(-2.3095135) q[0];
sx q[0];
rz(0.59659514) q[0];
rz(-2.3157628) q[1];
sx q[1];
rz(-1.700371) q[1];
sx q[1];
rz(1.9155496) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7694089) q[0];
sx q[0];
rz(-0.49159494) q[0];
sx q[0];
rz(0.43222897) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82310279) q[2];
sx q[2];
rz(-2.1777993) q[2];
sx q[2];
rz(-2.5603103) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.021124161) q[1];
sx q[1];
rz(-1.3509343) q[1];
sx q[1];
rz(0.0004417858) q[1];
rz(-pi) q[2];
rz(-0.0094718178) q[3];
sx q[3];
rz(-0.99820271) q[3];
sx q[3];
rz(-0.026281683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3423959) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(0.33102316) q[2];
rz(-2.3349169) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(1.4276918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1598635) q[0];
sx q[0];
rz(-1.9112497) q[0];
sx q[0];
rz(1.8925517) q[0];
rz(-3.0535835) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(2.1121315) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91988504) q[0];
sx q[0];
rz(-1.3516597) q[0];
sx q[0];
rz(-2.1973781) q[0];
rz(-2.7715893) q[2];
sx q[2];
rz(-0.68379935) q[2];
sx q[2];
rz(1.5497269) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0440812) q[1];
sx q[1];
rz(-1.71002) q[1];
sx q[1];
rz(-0.15686762) q[1];
rz(1.6920407) q[3];
sx q[3];
rz(-2.3070934) q[3];
sx q[3];
rz(0.27437011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.007894667) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(-0.50764817) q[2];
rz(1.7525904) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(-1.1631789) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4841109) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(1.5959651) q[0];
rz(-1.0428628) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(1.5159336) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7750054) q[0];
sx q[0];
rz(-0.89582743) q[0];
sx q[0];
rz(-1.7540027) q[0];
rz(-pi) q[1];
rz(3.0558673) q[2];
sx q[2];
rz(-2.1282196) q[2];
sx q[2];
rz(2.2103708) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.013853156) q[1];
sx q[1];
rz(-1.2989559) q[1];
sx q[1];
rz(-2.4747162) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58873119) q[3];
sx q[3];
rz(-1.7541459) q[3];
sx q[3];
rz(-2.6453032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3778014) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(-1.2949004) q[2];
rz(-3.0002248) q[3];
sx q[3];
rz(-0.54261345) q[3];
sx q[3];
rz(-0.98658371) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35448733) q[0];
sx q[0];
rz(-1.1802477) q[0];
sx q[0];
rz(1.6217344) q[0];
rz(0.63201085) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(2.246726) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8653523) q[0];
sx q[0];
rz(-1.9812752) q[0];
sx q[0];
rz(0.68666896) q[0];
rz(2.1331482) q[2];
sx q[2];
rz(-0.92698669) q[2];
sx q[2];
rz(1.2959727) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.289031) q[1];
sx q[1];
rz(-0.82419306) q[1];
sx q[1];
rz(-0.95552163) q[1];
x q[2];
rz(0.28646333) q[3];
sx q[3];
rz(-1.5781919) q[3];
sx q[3];
rz(-1.3163819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.52508369) q[2];
sx q[2];
rz(-2.775165) q[2];
sx q[2];
rz(1.1425225) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(-2.0660627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.557945) q[0];
sx q[0];
rz(-2.8201411) q[0];
sx q[0];
rz(-1.7640132) q[0];
rz(-0.47239834) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(-2.6766052) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87100055) q[0];
sx q[0];
rz(-2.5809079) q[0];
sx q[0];
rz(2.2339348) q[0];
rz(-pi) q[1];
rz(2.9372413) q[2];
sx q[2];
rz(-2.7933279) q[2];
sx q[2];
rz(-1.5262926) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7558414) q[1];
sx q[1];
rz(-2.5045536) q[1];
sx q[1];
rz(1.6511276) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70907866) q[3];
sx q[3];
rz(-1.7311829) q[3];
sx q[3];
rz(-1.3080314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.650699) q[2];
sx q[2];
rz(-1.8785672) q[2];
sx q[2];
rz(-0.4450376) q[2];
rz(2.2079091) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(0.26708189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0141107) q[0];
sx q[0];
rz(-1.575379) q[0];
sx q[0];
rz(1.7215464) q[0];
rz(3.1177915) q[1];
sx q[1];
rz(-2.5271466) q[1];
sx q[1];
rz(2.9856317) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1056846) q[0];
sx q[0];
rz(-1.3996291) q[0];
sx q[0];
rz(-0.25733421) q[0];
rz(-2.6290226) q[2];
sx q[2];
rz(-1.1935496) q[2];
sx q[2];
rz(-0.23114983) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5739463) q[1];
sx q[1];
rz(-1.5168377) q[1];
sx q[1];
rz(1.3957363) q[1];
x q[2];
rz(-2.6050623) q[3];
sx q[3];
rz(-2.0091669) q[3];
sx q[3];
rz(-2.8560672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0722787) q[2];
sx q[2];
rz(-0.57702714) q[2];
sx q[2];
rz(2.0689266) q[2];
rz(0.32564751) q[3];
sx q[3];
rz(-1.1762534) q[3];
sx q[3];
rz(1.5163039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1919365) q[0];
sx q[0];
rz(-0.40238109) q[0];
sx q[0];
rz(-2.8038213) q[0];
rz(2.0514964) q[1];
sx q[1];
rz(-2.1665159) q[1];
sx q[1];
rz(-2.8930194) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0524806) q[0];
sx q[0];
rz(-2.5626474) q[0];
sx q[0];
rz(0.46778932) q[0];
x q[1];
rz(-2.6612501) q[2];
sx q[2];
rz(-2.0311653) q[2];
sx q[2];
rz(-2.4127221) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.12411815) q[1];
sx q[1];
rz(-1.0867456) q[1];
sx q[1];
rz(-2.2494621) q[1];
rz(-pi) q[2];
rz(0.076898889) q[3];
sx q[3];
rz(-1.4014763) q[3];
sx q[3];
rz(2.3016735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8481855) q[2];
sx q[2];
rz(-2.6082787) q[2];
sx q[2];
rz(-0.08671134) q[2];
rz(-2.6596206) q[3];
sx q[3];
rz(-2.635699) q[3];
sx q[3];
rz(1.1530676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6329353) q[0];
sx q[0];
rz(-2.1338699) q[0];
sx q[0];
rz(-2.8588262) q[0];
rz(-0.70156082) q[1];
sx q[1];
rz(-0.82156721) q[1];
sx q[1];
rz(-1.823002) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1222262) q[0];
sx q[0];
rz(-1.1302233) q[0];
sx q[0];
rz(-3.0526572) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4852072) q[2];
sx q[2];
rz(-1.0644541) q[2];
sx q[2];
rz(-2.2040747) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1773771) q[1];
sx q[1];
rz(-2.4417158) q[1];
sx q[1];
rz(-1.3436951) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41140822) q[3];
sx q[3];
rz(-0.79380006) q[3];
sx q[3];
rz(-1.0994764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7302154) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(-0.39917699) q[2];
rz(-2.2579851) q[3];
sx q[3];
rz(-1.4727605) q[3];
sx q[3];
rz(-1.9201027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3783962) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(1.5989074) q[0];
rz(2.0762766) q[1];
sx q[1];
rz(-2.1677446) q[1];
sx q[1];
rz(1.261196) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3189145) q[0];
sx q[0];
rz(-2.7180674) q[0];
sx q[0];
rz(-0.25869297) q[0];
x q[1];
rz(0.59826675) q[2];
sx q[2];
rz(-1.5891799) q[2];
sx q[2];
rz(-2.604904) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8775455) q[1];
sx q[1];
rz(-2.0685158) q[1];
sx q[1];
rz(1.643814) q[1];
x q[2];
rz(-0.68393771) q[3];
sx q[3];
rz(-1.4468907) q[3];
sx q[3];
rz(-0.77139664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5836872) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(-2.5718001) q[2];
rz(-1.2184881) q[3];
sx q[3];
rz(-0.64703882) q[3];
sx q[3];
rz(-0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31310836) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(-0.60824153) q[1];
sx q[1];
rz(-0.47641644) q[1];
sx q[1];
rz(0.48412916) q[1];
rz(2.1424978) q[2];
sx q[2];
rz(-1.5487164) q[2];
sx q[2];
rz(-1.6185417) q[2];
rz(-1.7096056) q[3];
sx q[3];
rz(-2.1574253) q[3];
sx q[3];
rz(0.16711259) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
