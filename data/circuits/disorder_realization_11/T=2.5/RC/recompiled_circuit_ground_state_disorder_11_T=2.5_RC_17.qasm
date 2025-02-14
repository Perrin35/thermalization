OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0271072) q[0];
sx q[0];
rz(-3.0740102) q[0];
sx q[0];
rz(-2.5525868) q[0];
rz(2.0677805) q[1];
sx q[1];
rz(-1.685073) q[1];
sx q[1];
rz(-2.9340802) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1197101) q[0];
sx q[0];
rz(-1.8069022) q[0];
sx q[0];
rz(-2.9364999) q[0];
rz(-pi) q[1];
rz(1.6512538) q[2];
sx q[2];
rz(-1.7535889) q[2];
sx q[2];
rz(-1.8535943) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.0879651) q[1];
sx q[1];
rz(-0.028465406) q[1];
sx q[1];
rz(1.444412) q[1];
x q[2];
rz(-0.46300827) q[3];
sx q[3];
rz(-1.6673458) q[3];
sx q[3];
rz(-2.2322886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33240685) q[2];
sx q[2];
rz(-0.016533479) q[2];
sx q[2];
rz(0.22745505) q[2];
rz(-2.4849232) q[3];
sx q[3];
rz(-2.4469817) q[3];
sx q[3];
rz(-2.606126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4559795) q[0];
sx q[0];
rz(-3.1094636) q[0];
sx q[0];
rz(2.6973714) q[0];
rz(-1.1086858) q[1];
sx q[1];
rz(-1.8007092) q[1];
sx q[1];
rz(1.1681555) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.927181) q[0];
sx q[0];
rz(-2.4200038) q[0];
sx q[0];
rz(2.847228) q[0];
rz(-pi) q[1];
rz(0.45380317) q[2];
sx q[2];
rz(-1.6112932) q[2];
sx q[2];
rz(-1.7292505) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.44732783) q[1];
sx q[1];
rz(-0.45644293) q[1];
sx q[1];
rz(1.3988628) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1862964) q[3];
sx q[3];
rz(-1.9855193) q[3];
sx q[3];
rz(0.44594582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4593959) q[2];
sx q[2];
rz(-1.0707868) q[2];
sx q[2];
rz(0.11967858) q[2];
rz(-2.1648572) q[3];
sx q[3];
rz(-0.027438199) q[3];
sx q[3];
rz(-1.1915092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6617436) q[0];
sx q[0];
rz(-2.9763344) q[0];
sx q[0];
rz(-1.6462434) q[0];
rz(-2.7961075) q[1];
sx q[1];
rz(-0.59436878) q[1];
sx q[1];
rz(-0.77846175) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0986106) q[0];
sx q[0];
rz(-1.6263282) q[0];
sx q[0];
rz(1.5639202) q[0];
x q[1];
rz(2.5359306) q[2];
sx q[2];
rz(-2.2251943) q[2];
sx q[2];
rz(-0.6336279) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.045588551) q[1];
sx q[1];
rz(-2.0109777) q[1];
sx q[1];
rz(0.11543302) q[1];
rz(2.4308079) q[3];
sx q[3];
rz(-0.87276283) q[3];
sx q[3];
rz(-0.16434114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8567132) q[2];
sx q[2];
rz(-3.1046125) q[2];
sx q[2];
rz(1.7516837) q[2];
rz(-3.1189647) q[3];
sx q[3];
rz(-0.32630625) q[3];
sx q[3];
rz(2.7271395) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49700272) q[0];
sx q[0];
rz(-3.1019326) q[0];
sx q[0];
rz(0.52763754) q[0];
rz(-0.35609326) q[1];
sx q[1];
rz(-1.6368607) q[1];
sx q[1];
rz(1.5632695) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.281412) q[0];
sx q[0];
rz(-1.9989655) q[0];
sx q[0];
rz(-0.58089818) q[0];
rz(-3.061297) q[2];
sx q[2];
rz(-0.67932898) q[2];
sx q[2];
rz(1.7905362) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4112579) q[1];
sx q[1];
rz(-2.6739542) q[1];
sx q[1];
rz(-1.1958666) q[1];
rz(-pi) q[2];
rz(-2.0464055) q[3];
sx q[3];
rz(-1.4970253) q[3];
sx q[3];
rz(0.99223415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.57033527) q[2];
sx q[2];
rz(-3.1331077) q[2];
sx q[2];
rz(-0.045489475) q[2];
rz(2.3174543) q[3];
sx q[3];
rz(-2.6233311) q[3];
sx q[3];
rz(0.99388188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5387886) q[0];
sx q[0];
rz(-2.0863918) q[0];
sx q[0];
rz(-1.3987199) q[0];
rz(2.973373) q[1];
sx q[1];
rz(-1.8054447) q[1];
sx q[1];
rz(2.9199563) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0960064) q[0];
sx q[0];
rz(-1.0500589) q[0];
sx q[0];
rz(-0.27294275) q[0];
rz(-pi) q[1];
rz(3.0732642) q[2];
sx q[2];
rz(-1.7842275) q[2];
sx q[2];
rz(-2.1377856) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6252215) q[1];
sx q[1];
rz(-2.0091426) q[1];
sx q[1];
rz(-0.0268374) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6898674) q[3];
sx q[3];
rz(-2.171031) q[3];
sx q[3];
rz(-0.33590318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98442709) q[2];
sx q[2];
rz(-0.027313622) q[2];
sx q[2];
rz(-2.1402764) q[2];
rz(2.2973513) q[3];
sx q[3];
rz(-0.068035754) q[3];
sx q[3];
rz(2.3943353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.35185128) q[0];
sx q[0];
rz(-0.25775596) q[0];
sx q[0];
rz(-1.7716273) q[0];
rz(-0.17579707) q[1];
sx q[1];
rz(-1.628123) q[1];
sx q[1];
rz(1.0584077) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8254614) q[0];
sx q[0];
rz(-0.93535813) q[0];
sx q[0];
rz(-2.7881289) q[0];
rz(0.01137017) q[2];
sx q[2];
rz(-0.98207131) q[2];
sx q[2];
rz(-2.7259072) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0479916) q[1];
sx q[1];
rz(-2.2696884) q[1];
sx q[1];
rz(1.2456139) q[1];
rz(-pi) q[2];
rz(-2.2203683) q[3];
sx q[3];
rz(-1.2787673) q[3];
sx q[3];
rz(2.667922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.83219641) q[2];
sx q[2];
rz(-3.1377073) q[2];
sx q[2];
rz(-0.84581172) q[2];
rz(-2.5064365) q[3];
sx q[3];
rz(-2.784415) q[3];
sx q[3];
rz(0.14491189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.7307067) q[0];
sx q[0];
rz(-2.9697953) q[0];
sx q[0];
rz(0.26096499) q[0];
rz(-1.4313401) q[1];
sx q[1];
rz(-0.15161082) q[1];
sx q[1];
rz(1.8401015) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8943202) q[0];
sx q[0];
rz(-1.1218583) q[0];
sx q[0];
rz(-2.2821064) q[0];
rz(-pi) q[1];
rz(-0.11110961) q[2];
sx q[2];
rz(-1.5556364) q[2];
sx q[2];
rz(0.76778256) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.83205081) q[1];
sx q[1];
rz(-2.0027035) q[1];
sx q[1];
rz(-1.5276272) q[1];
rz(-1.8347307) q[3];
sx q[3];
rz(-1.5698877) q[3];
sx q[3];
rz(3.0543882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5362376) q[2];
sx q[2];
rz(-1.321512) q[2];
sx q[2];
rz(3.1129254) q[2];
rz(0.043449314) q[3];
sx q[3];
rz(-2.9048007) q[3];
sx q[3];
rz(-1.4437599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0726149) q[0];
sx q[0];
rz(-0.34729877) q[0];
sx q[0];
rz(2.8436227) q[0];
rz(-1.7295674) q[1];
sx q[1];
rz(-2.7822918) q[1];
sx q[1];
rz(-1.5922458) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6291154) q[0];
sx q[0];
rz(-0.40226679) q[0];
sx q[0];
rz(0.8921807) q[0];
rz(-3.1352714) q[2];
sx q[2];
rz(-1.5875193) q[2];
sx q[2];
rz(-1.4249546) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(5/(7*pi)) q[1];
sx q[1];
rz(-0.67247219) q[1];
sx q[1];
rz(-0.48788957) q[1];
rz(-pi) q[2];
rz(-1.8531591) q[3];
sx q[3];
rz(-1.3579277) q[3];
sx q[3];
rz(-2.4119056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2397502) q[2];
sx q[2];
rz(-0.032278927) q[2];
sx q[2];
rz(0.65995222) q[2];
rz(0.7800855) q[3];
sx q[3];
rz(-1.2939021) q[3];
sx q[3];
rz(1.1680781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50769794) q[0];
sx q[0];
rz(-0.036490353) q[0];
sx q[0];
rz(-2.3376035) q[0];
rz(-2.9102303) q[1];
sx q[1];
rz(-1.7295126) q[1];
sx q[1];
rz(0.075459935) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20310878) q[0];
sx q[0];
rz(-1.9265947) q[0];
sx q[0];
rz(1.7901161) q[0];
x q[1];
rz(3.6663975e-05) q[2];
sx q[2];
rz(-1.5701887) q[2];
sx q[2];
rz(-3.031784) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.67565221) q[1];
sx q[1];
rz(-2.2237101) q[1];
sx q[1];
rz(0.5799579) q[1];
rz(0.034377957) q[3];
sx q[3];
rz(-1.9581736) q[3];
sx q[3];
rz(2.9971214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1239473) q[2];
sx q[2];
rz(-0.46111527) q[2];
sx q[2];
rz(-0.41575113) q[2];
rz(1.4990643) q[3];
sx q[3];
rz(-3.1406904) q[3];
sx q[3];
rz(-0.7846964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42522624) q[0];
sx q[0];
rz(-3.0854736) q[0];
sx q[0];
rz(2.8477493) q[0];
rz(1.6055239) q[1];
sx q[1];
rz(-0.87296456) q[1];
sx q[1];
rz(1.4651089) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11742453) q[0];
sx q[0];
rz(-3.0367622) q[0];
sx q[0];
rz(0.094314055) q[0];
rz(2.5757786) q[2];
sx q[2];
rz(-2.4560258) q[2];
sx q[2];
rz(-3.1367347) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5398354) q[1];
sx q[1];
rz(-2.4944427) q[1];
sx q[1];
rz(-0.045250968) q[1];
x q[2];
rz(2.5876643) q[3];
sx q[3];
rz(-3.1331473) q[3];
sx q[3];
rz(2.3609207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6426223) q[2];
sx q[2];
rz(-0.62752807) q[2];
sx q[2];
rz(1.2205396) q[2];
rz(-2.7484861) q[3];
sx q[3];
rz(-0.016963907) q[3];
sx q[3];
rz(-0.19256798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6846234) q[0];
sx q[0];
rz(-1.715281) q[0];
sx q[0];
rz(2.793978) q[0];
rz(0.63812232) q[1];
sx q[1];
rz(-0.77824021) q[1];
sx q[1];
rz(2.4660769) q[1];
rz(-0.19370041) q[2];
sx q[2];
rz(-1.4932409) q[2];
sx q[2];
rz(1.3328339) q[2];
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
