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
rz(-1.00296) q[0];
sx q[0];
rz(-2.668219) q[0];
sx q[0];
rz(1.0869429) q[0];
rz(2.3776157) q[1];
sx q[1];
rz(-2.6978701) q[1];
sx q[1];
rz(0.8313764) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2646541) q[0];
sx q[0];
rz(-2.0994517) q[0];
sx q[0];
rz(-1.009788) q[0];
x q[1];
rz(-0.23687266) q[2];
sx q[2];
rz(-2.028596) q[2];
sx q[2];
rz(-1.2352236) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9688111) q[1];
sx q[1];
rz(-2.9829142) q[1];
sx q[1];
rz(-2.2814343) q[1];
rz(-pi) q[2];
x q[2];
rz(0.071962996) q[3];
sx q[3];
rz(-1.5974034) q[3];
sx q[3];
rz(-3.0940477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6143371) q[2];
sx q[2];
rz(-1.273512) q[2];
sx q[2];
rz(-2.7296076) q[2];
rz(-2.8968503) q[3];
sx q[3];
rz(-2.5089846) q[3];
sx q[3];
rz(0.28732029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43993846) q[0];
sx q[0];
rz(-2.169862) q[0];
sx q[0];
rz(0.42020759) q[0];
rz(2.6523051) q[1];
sx q[1];
rz(-2.2261765) q[1];
sx q[1];
rz(-1.9583826) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1294043) q[0];
sx q[0];
rz(-2.3478635) q[0];
sx q[0];
rz(-0.94214006) q[0];
rz(-0.80090307) q[2];
sx q[2];
rz(-2.704853) q[2];
sx q[2];
rz(0.85608053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4192761) q[1];
sx q[1];
rz(-2.4240383) q[1];
sx q[1];
rz(0.60075398) q[1];
rz(-1.39721) q[3];
sx q[3];
rz(-1.6866444) q[3];
sx q[3];
rz(-0.66044745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.45541397) q[2];
sx q[2];
rz(-1.4718461) q[2];
sx q[2];
rz(0.14032826) q[2];
rz(-0.27648196) q[3];
sx q[3];
rz(-0.77767196) q[3];
sx q[3];
rz(-2.8592143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33271933) q[0];
sx q[0];
rz(-2.7820899) q[0];
sx q[0];
rz(-1.7669539) q[0];
rz(-0.91284347) q[1];
sx q[1];
rz(-2.5990867) q[1];
sx q[1];
rz(1.0309781) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1682273) q[0];
sx q[0];
rz(-1.7534257) q[0];
sx q[0];
rz(-1.2398861) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1314083) q[2];
sx q[2];
rz(-1.2282145) q[2];
sx q[2];
rz(-0.48784384) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.678434) q[1];
sx q[1];
rz(-2.7292876) q[1];
sx q[1];
rz(-1.5458471) q[1];
rz(-2.8171478) q[3];
sx q[3];
rz(-1.3881171) q[3];
sx q[3];
rz(2.7871391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1294407) q[2];
sx q[2];
rz(-2.0871711) q[2];
sx q[2];
rz(-0.63808179) q[2];
rz(-0.10073999) q[3];
sx q[3];
rz(-2.1295857) q[3];
sx q[3];
rz(0.49873275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8345555) q[0];
sx q[0];
rz(-1.6224253) q[0];
sx q[0];
rz(1.3822973) q[0];
rz(-2.8221829) q[1];
sx q[1];
rz(-1.5089792) q[1];
sx q[1];
rz(0.75739783) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46947458) q[0];
sx q[0];
rz(-1.645477) q[0];
sx q[0];
rz(-1.5622219) q[0];
x q[1];
rz(1.6998197) q[2];
sx q[2];
rz(-1.4296544) q[2];
sx q[2];
rz(-0.8069969) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3925303) q[1];
sx q[1];
rz(-0.99732256) q[1];
sx q[1];
rz(-0.5368781) q[1];
rz(-2.4457127) q[3];
sx q[3];
rz(-0.73344389) q[3];
sx q[3];
rz(-0.91737285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.77183977) q[2];
sx q[2];
rz(-2.2660793) q[2];
sx q[2];
rz(0.80779752) q[2];
rz(1.7317023) q[3];
sx q[3];
rz(-1.2597224) q[3];
sx q[3];
rz(0.65214777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0659502) q[0];
sx q[0];
rz(-0.37794161) q[0];
sx q[0];
rz(-2.156303) q[0];
rz(1.6458884) q[1];
sx q[1];
rz(-2.0031877) q[1];
sx q[1];
rz(2.2122502) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2332704) q[0];
sx q[0];
rz(-2.4886971) q[0];
sx q[0];
rz(-2.6920431) q[0];
rz(-pi) q[1];
rz(1.1628289) q[2];
sx q[2];
rz(-1.0024602) q[2];
sx q[2];
rz(1.6108244) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.73354355) q[1];
sx q[1];
rz(-1.0449303) q[1];
sx q[1];
rz(-2.9108564) q[1];
rz(-pi) q[2];
rz(2.1319904) q[3];
sx q[3];
rz(-1.5586886) q[3];
sx q[3];
rz(0.59134134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.052496584) q[2];
sx q[2];
rz(-1.7868944) q[2];
sx q[2];
rz(2.055114) q[2];
rz(2.1975482) q[3];
sx q[3];
rz(-3.0510674) q[3];
sx q[3];
rz(-2.0199147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97359598) q[0];
sx q[0];
rz(-2.7041628) q[0];
sx q[0];
rz(2.9214389) q[0];
rz(1.6379697) q[1];
sx q[1];
rz(-1.3238246) q[1];
sx q[1];
rz(-0.55353177) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39237994) q[0];
sx q[0];
rz(-0.22870453) q[0];
sx q[0];
rz(-2.6116651) q[0];
x q[1];
rz(-2.7017864) q[2];
sx q[2];
rz(-0.39798073) q[2];
sx q[2];
rz(0.27962886) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.625861) q[1];
sx q[1];
rz(-1.5222094) q[1];
sx q[1];
rz(-0.69339417) q[1];
x q[2];
rz(-1.342085) q[3];
sx q[3];
rz(-0.89806496) q[3];
sx q[3];
rz(-1.7941689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0963875) q[2];
sx q[2];
rz(-1.4676981) q[2];
sx q[2];
rz(-0.20827797) q[2];
rz(-2.0221209) q[3];
sx q[3];
rz(-1.2866311) q[3];
sx q[3];
rz(-2.313405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.646362) q[0];
sx q[0];
rz(-1.3883075) q[0];
sx q[0];
rz(2.0679423) q[0];
rz(3.011009) q[1];
sx q[1];
rz(-1.5675631) q[1];
sx q[1];
rz(-1.0098339) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9511418) q[0];
sx q[0];
rz(-1.9660239) q[0];
sx q[0];
rz(1.7717096) q[0];
rz(-pi) q[1];
rz(-2.9931941) q[2];
sx q[2];
rz(-1.6434533) q[2];
sx q[2];
rz(2.1138482) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1853789) q[1];
sx q[1];
rz(-1.662743) q[1];
sx q[1];
rz(-2.5673709) q[1];
rz(-pi) q[2];
rz(0.66422446) q[3];
sx q[3];
rz(-1.281732) q[3];
sx q[3];
rz(-2.5632896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4449571) q[2];
sx q[2];
rz(-0.37142763) q[2];
sx q[2];
rz(-2.8182287) q[2];
rz(0.67445406) q[3];
sx q[3];
rz(-0.89265299) q[3];
sx q[3];
rz(-3.118012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12482878) q[0];
sx q[0];
rz(-2.6478196) q[0];
sx q[0];
rz(1.0892518) q[0];
rz(-0.59843868) q[1];
sx q[1];
rz(-1.2804223) q[1];
sx q[1];
rz(-2.3425897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.625811) q[0];
sx q[0];
rz(-2.7220753) q[0];
sx q[0];
rz(-2.3965492) q[0];
rz(-pi) q[1];
rz(2.2369821) q[2];
sx q[2];
rz(-1.4965726) q[2];
sx q[2];
rz(0.42236172) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.028513718) q[1];
sx q[1];
rz(-1.4220107) q[1];
sx q[1];
rz(2.8362175) q[1];
rz(-pi) q[2];
x q[2];
rz(1.576242) q[3];
sx q[3];
rz(-1.8126235) q[3];
sx q[3];
rz(-0.33978811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4716855) q[2];
sx q[2];
rz(-0.50614637) q[2];
sx q[2];
rz(-1.252582) q[2];
rz(1.0505098) q[3];
sx q[3];
rz(-2.551008) q[3];
sx q[3];
rz(-2.2091776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2893386) q[0];
sx q[0];
rz(-1.4006571) q[0];
sx q[0];
rz(0.5156714) q[0];
rz(-1.4562666) q[1];
sx q[1];
rz(-1.8003502) q[1];
sx q[1];
rz(-2.8175443) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0312248) q[0];
sx q[0];
rz(-0.47981167) q[0];
sx q[0];
rz(-2.5672002) q[0];
x q[1];
rz(-2.134974) q[2];
sx q[2];
rz(-0.7340275) q[2];
sx q[2];
rz(-0.013040868) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9726213) q[1];
sx q[1];
rz(-1.4310799) q[1];
sx q[1];
rz(-2.0129544) q[1];
x q[2];
rz(1.4435931) q[3];
sx q[3];
rz(-2.5802744) q[3];
sx q[3];
rz(1.4385245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1590165) q[2];
sx q[2];
rz(-1.1218718) q[2];
sx q[2];
rz(-2.4328361) q[2];
rz(2.3580264) q[3];
sx q[3];
rz(-1.9400027) q[3];
sx q[3];
rz(1.5172575) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5358955) q[0];
sx q[0];
rz(-0.63740969) q[0];
sx q[0];
rz(-2.9851483) q[0];
rz(-2.5018196) q[1];
sx q[1];
rz(-2.4867058) q[1];
sx q[1];
rz(-2.1515062) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7759686) q[0];
sx q[0];
rz(-2.1556615) q[0];
sx q[0];
rz(0.82230277) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9205356) q[2];
sx q[2];
rz(-1.4863401) q[2];
sx q[2];
rz(-1.3049558) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1875547) q[1];
sx q[1];
rz(-2.1219664) q[1];
sx q[1];
rz(-1.7269082) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5699432) q[3];
sx q[3];
rz(-0.23861966) q[3];
sx q[3];
rz(2.1941136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5797552) q[2];
sx q[2];
rz(-1.729915) q[2];
sx q[2];
rz(-2.5282395) q[2];
rz(0.26398811) q[3];
sx q[3];
rz(-1.3365021) q[3];
sx q[3];
rz(-0.67155513) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092125208) q[0];
sx q[0];
rz(-1.5277852) q[0];
sx q[0];
rz(1.1396136) q[0];
rz(-1.4641948) q[1];
sx q[1];
rz(-0.72036998) q[1];
sx q[1];
rz(1.5710685) q[1];
rz(1.9665267) q[2];
sx q[2];
rz(-1.1968975) q[2];
sx q[2];
rz(0.71002985) q[2];
rz(3.0861978) q[3];
sx q[3];
rz(-2.5221586) q[3];
sx q[3];
rz(0.4172162) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
