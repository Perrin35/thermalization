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
rz(-1.3396076) q[0];
sx q[0];
rz(2.7924502) q[0];
sx q[0];
rz(10.452441) q[0];
rz(-0.034962058) q[1];
sx q[1];
rz(-1.0171913) q[1];
sx q[1];
rz(-1.4453759) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48533121) q[0];
sx q[0];
rz(-0.6888656) q[0];
sx q[0];
rz(0.553225) q[0];
rz(-1.6048172) q[2];
sx q[2];
rz(-1.3933225) q[2];
sx q[2];
rz(-0.79171514) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.45856341) q[1];
sx q[1];
rz(-1.7818319) q[1];
sx q[1];
rz(0.33331916) q[1];
rz(2.7463116) q[3];
sx q[3];
rz(-2.3594405) q[3];
sx q[3];
rz(-0.30667337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8029636) q[2];
sx q[2];
rz(-2.6821319) q[2];
sx q[2];
rz(0.53172338) q[2];
rz(-1.3945329) q[3];
sx q[3];
rz(-1.2732384) q[3];
sx q[3];
rz(1.9664221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26674536) q[0];
sx q[0];
rz(-1.6371472) q[0];
sx q[0];
rz(-2.3496085) q[0];
rz(-2.2199471) q[1];
sx q[1];
rz(-1.4896723) q[1];
sx q[1];
rz(2.0335061) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32852325) q[0];
sx q[0];
rz(-1.6824772) q[0];
sx q[0];
rz(-0.96891667) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94723786) q[2];
sx q[2];
rz(-1.6102566) q[2];
sx q[2];
rz(0.90496162) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0597569) q[1];
sx q[1];
rz(-1.949777) q[1];
sx q[1];
rz(2.3339218) q[1];
rz(1.2543801) q[3];
sx q[3];
rz(-1.5174447) q[3];
sx q[3];
rz(-0.90094588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.38829982) q[2];
sx q[2];
rz(-2.1123977) q[2];
sx q[2];
rz(1.252906) q[2];
rz(0.62475723) q[3];
sx q[3];
rz(-1.2478991) q[3];
sx q[3];
rz(-2.6784082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4598292) q[0];
sx q[0];
rz(-2.9637931) q[0];
sx q[0];
rz(-2.8681927) q[0];
rz(3.0699442) q[1];
sx q[1];
rz(-1.1337846) q[1];
sx q[1];
rz(1.515548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1119627) q[0];
sx q[0];
rz(-0.66034895) q[0];
sx q[0];
rz(2.4899954) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63457527) q[2];
sx q[2];
rz(-2.5800642) q[2];
sx q[2];
rz(0.47088366) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7715986) q[1];
sx q[1];
rz(-0.99111667) q[1];
sx q[1];
rz(1.3808151) q[1];
x q[2];
rz(-2.1642926) q[3];
sx q[3];
rz(-1.7827099) q[3];
sx q[3];
rz(-2.7023072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7800954) q[2];
sx q[2];
rz(-0.7146892) q[2];
sx q[2];
rz(1.3698461) q[2];
rz(-2.0731549) q[3];
sx q[3];
rz(-1.7504102) q[3];
sx q[3];
rz(0.95503241) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84900981) q[0];
sx q[0];
rz(-2.7136901) q[0];
sx q[0];
rz(-3.0191315) q[0];
rz(0.017008688) q[1];
sx q[1];
rz(-1.5127134) q[1];
sx q[1];
rz(0.88517991) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6017319) q[0];
sx q[0];
rz(-2.0575876) q[0];
sx q[0];
rz(-0.15386015) q[0];
rz(-pi) q[1];
rz(0.83577581) q[2];
sx q[2];
rz(-1.3124497) q[2];
sx q[2];
rz(1.5193743) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0065476848) q[1];
sx q[1];
rz(-1.1741877) q[1];
sx q[1];
rz(2.066925) q[1];
rz(1.7649129) q[3];
sx q[3];
rz(-2.2633865) q[3];
sx q[3];
rz(3.0339981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7077606) q[2];
sx q[2];
rz(-1.3154575) q[2];
sx q[2];
rz(-3.0641277) q[2];
rz(-0.42099434) q[3];
sx q[3];
rz(-1.1869895) q[3];
sx q[3];
rz(-0.25119701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0780599) q[0];
sx q[0];
rz(-0.01883004) q[0];
sx q[0];
rz(-0.68191648) q[0];
rz(-0.49184999) q[1];
sx q[1];
rz(-1.0268772) q[1];
sx q[1];
rz(-1.9815725) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9043961) q[0];
sx q[0];
rz(-0.86441308) q[0];
sx q[0];
rz(-2.893413) q[0];
x q[1];
rz(-3.0536041) q[2];
sx q[2];
rz(-2.5694429) q[2];
sx q[2];
rz(2.3038626) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.3403389) q[1];
sx q[1];
rz(-1.7427708) q[1];
sx q[1];
rz(-0.94229631) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8468708) q[3];
sx q[3];
rz(-1.1663166) q[3];
sx q[3];
rz(-1.3724788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.68449768) q[2];
sx q[2];
rz(-0.86898494) q[2];
sx q[2];
rz(2.1333466) q[2];
rz(-1.5021987) q[3];
sx q[3];
rz(-2.0366171) q[3];
sx q[3];
rz(-0.26386279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0428001) q[0];
sx q[0];
rz(-1.5473939) q[0];
sx q[0];
rz(0.40476009) q[0];
rz(-0.81819355) q[1];
sx q[1];
rz(-1.300756) q[1];
sx q[1];
rz(2.4249446) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34490624) q[0];
sx q[0];
rz(-2.8622238) q[0];
sx q[0];
rz(-3.1059103) q[0];
rz(-pi) q[1];
rz(-1.0425074) q[2];
sx q[2];
rz(-2.2620438) q[2];
sx q[2];
rz(1.3529568) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.055775) q[1];
sx q[1];
rz(-1.1113864) q[1];
sx q[1];
rz(1.8618078) q[1];
rz(-1.9394373) q[3];
sx q[3];
rz(-1.0615942) q[3];
sx q[3];
rz(-0.64083767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8019668) q[2];
sx q[2];
rz(-0.64029396) q[2];
sx q[2];
rz(-0.65529811) q[2];
rz(0.35704923) q[3];
sx q[3];
rz(-1.3909631) q[3];
sx q[3];
rz(-2.6220139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.201467) q[0];
sx q[0];
rz(-2.8253912) q[0];
sx q[0];
rz(-1.0200208) q[0];
rz(0.54234281) q[1];
sx q[1];
rz(-2.0261363) q[1];
sx q[1];
rz(-1.382359) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6868149) q[0];
sx q[0];
rz(-2.1612148) q[0];
sx q[0];
rz(-0.83417828) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7075536) q[2];
sx q[2];
rz(-1.7801577) q[2];
sx q[2];
rz(-0.48780045) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2445557) q[1];
sx q[1];
rz(-0.4677597) q[1];
sx q[1];
rz(-1.7451343) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1145688) q[3];
sx q[3];
rz(-3.0579429) q[3];
sx q[3];
rz(0.25318709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0688087) q[2];
sx q[2];
rz(-1.5280318) q[2];
sx q[2];
rz(-2.7911348) q[2];
rz(0.46640486) q[3];
sx q[3];
rz(-0.91752183) q[3];
sx q[3];
rz(-2.1980227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7826409) q[0];
sx q[0];
rz(-0.39059165) q[0];
sx q[0];
rz(2.1851831) q[0];
rz(1.6920754) q[1];
sx q[1];
rz(-2.3608975) q[1];
sx q[1];
rz(-2.4704959) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22204493) q[0];
sx q[0];
rz(-1.4853108) q[0];
sx q[0];
rz(1.4287022) q[0];
rz(-pi) q[1];
rz(1.7809733) q[2];
sx q[2];
rz(-2.1297014) q[2];
sx q[2];
rz(2.0483537) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.33167142) q[1];
sx q[1];
rz(-1.5070032) q[1];
sx q[1];
rz(1.924519) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8160041) q[3];
sx q[3];
rz(-0.88174654) q[3];
sx q[3];
rz(-2.6264555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2676919) q[2];
sx q[2];
rz(-2.716422) q[2];
sx q[2];
rz(1.1262061) q[2];
rz(-2.8113484) q[3];
sx q[3];
rz(-1.4010022) q[3];
sx q[3];
rz(-0.11708524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0432334) q[0];
sx q[0];
rz(-0.57906228) q[0];
sx q[0];
rz(3.0845355) q[0];
rz(3.0077899) q[1];
sx q[1];
rz(-1.0452784) q[1];
sx q[1];
rz(0.71279508) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4445515) q[0];
sx q[0];
rz(-1.6336055) q[0];
sx q[0];
rz(2.8256326) q[0];
rz(-pi) q[1];
rz(1.1634689) q[2];
sx q[2];
rz(-1.1962657) q[2];
sx q[2];
rz(2.5782812) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9660898) q[1];
sx q[1];
rz(-2.0398629) q[1];
sx q[1];
rz(-1.1203946) q[1];
rz(-pi) q[2];
rz(1.6678048) q[3];
sx q[3];
rz(-0.69655124) q[3];
sx q[3];
rz(-2.4318621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.628525) q[2];
sx q[2];
rz(-1.2863938) q[2];
sx q[2];
rz(-0.081550278) q[2];
rz(1.2201803) q[3];
sx q[3];
rz(-2.8704075) q[3];
sx q[3];
rz(1.0903821) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9995025) q[0];
sx q[0];
rz(-1.5220078) q[0];
sx q[0];
rz(-1.5600486) q[0];
rz(-1.1355431) q[1];
sx q[1];
rz(-1.138843) q[1];
sx q[1];
rz(-1.9948237) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1449908) q[0];
sx q[0];
rz(-1.0787316) q[0];
sx q[0];
rz(-0.73623118) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3421971) q[2];
sx q[2];
rz(-2.172432) q[2];
sx q[2];
rz(2.5863199) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.15849491) q[1];
sx q[1];
rz(-1.4870411) q[1];
sx q[1];
rz(0.2106481) q[1];
rz(1.7079748) q[3];
sx q[3];
rz(-0.56932031) q[3];
sx q[3];
rz(-2.5018321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5494988) q[2];
sx q[2];
rz(-1.8213976) q[2];
sx q[2];
rz(2.7756694) q[2];
rz(-0.38771114) q[3];
sx q[3];
rz(-1.1185027) q[3];
sx q[3];
rz(1.6754735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.6112919) q[0];
sx q[0];
rz(-0.74321754) q[0];
sx q[0];
rz(2.8331953) q[0];
rz(0.74956924) q[1];
sx q[1];
rz(-2.0105965) q[1];
sx q[1];
rz(-1.2731193) q[1];
rz(-0.34229924) q[2];
sx q[2];
rz(-2.2020349) q[2];
sx q[2];
rz(2.0943691) q[2];
rz(-0.6324296) q[3];
sx q[3];
rz(-1.2865744) q[3];
sx q[3];
rz(-2.1435973) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
