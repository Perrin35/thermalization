OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.95076686) q[0];
sx q[0];
rz(-1.8523676) q[0];
sx q[0];
rz(0.092913203) q[0];
rz(-0.095280401) q[1];
sx q[1];
rz(-0.73268259) q[1];
sx q[1];
rz(1.9021775) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5372972) q[0];
sx q[0];
rz(-1.822246) q[0];
sx q[0];
rz(0.23668134) q[0];
x q[1];
rz(1.7669414) q[2];
sx q[2];
rz(-1.8298379) q[2];
sx q[2];
rz(2.1262622) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.8653633) q[1];
sx q[1];
rz(-1.6003834) q[1];
sx q[1];
rz(2.8113356) q[1];
rz(-2.5704424) q[3];
sx q[3];
rz(-1.1451654) q[3];
sx q[3];
rz(-1.7895123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6713082) q[2];
sx q[2];
rz(-1.4802063) q[2];
sx q[2];
rz(1.3489464) q[2];
rz(0.74364439) q[3];
sx q[3];
rz(-2.8642004) q[3];
sx q[3];
rz(-2.9644137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.141356) q[0];
sx q[0];
rz(-1.5445222) q[0];
sx q[0];
rz(0.29362383) q[0];
rz(1.0307505) q[1];
sx q[1];
rz(-1.7799957) q[1];
sx q[1];
rz(-0.34293276) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6688419) q[0];
sx q[0];
rz(-1.051051) q[0];
sx q[0];
rz(0.23730554) q[0];
x q[1];
rz(0.12983506) q[2];
sx q[2];
rz(-1.9399683) q[2];
sx q[2];
rz(2.2463727) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56655069) q[1];
sx q[1];
rz(-2.2600265) q[1];
sx q[1];
rz(-1.8495967) q[1];
x q[2];
rz(-0.94128709) q[3];
sx q[3];
rz(-1.3664748) q[3];
sx q[3];
rz(-2.2245537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.12884101) q[2];
sx q[2];
rz(-1.7495456) q[2];
sx q[2];
rz(-2.3823605) q[2];
rz(-3.1208842) q[3];
sx q[3];
rz(-1.2660675) q[3];
sx q[3];
rz(2.7032963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52221209) q[0];
sx q[0];
rz(-0.601957) q[0];
sx q[0];
rz(-3.1371327) q[0];
rz(-0.64747539) q[1];
sx q[1];
rz(-0.59440333) q[1];
sx q[1];
rz(2.3562145) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68092184) q[0];
sx q[0];
rz(-0.17853949) q[0];
sx q[0];
rz(0.4692678) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.17526971) q[2];
sx q[2];
rz(-1.2614377) q[2];
sx q[2];
rz(-0.64495211) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1238034) q[1];
sx q[1];
rz(-0.54517704) q[1];
sx q[1];
rz(1.4383609) q[1];
rz(2.3397869) q[3];
sx q[3];
rz(-2.2983645) q[3];
sx q[3];
rz(-1.8319195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7711827) q[2];
sx q[2];
rz(-2.2213171) q[2];
sx q[2];
rz(1.3578337) q[2];
rz(-0.34058288) q[3];
sx q[3];
rz(-0.95652306) q[3];
sx q[3];
rz(2.3497439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(-0.99712813) q[0];
sx q[0];
rz(-2.0576532) q[0];
sx q[0];
rz(2.2564364) q[0];
rz(-0.74686933) q[1];
sx q[1];
rz(-0.62364686) q[1];
sx q[1];
rz(-1.0964099) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80712705) q[0];
sx q[0];
rz(-1.8821446) q[0];
sx q[0];
rz(1.2047355) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6302413) q[2];
sx q[2];
rz(-1.5402093) q[2];
sx q[2];
rz(-0.93418834) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4064286) q[1];
sx q[1];
rz(-0.53710912) q[1];
sx q[1];
rz(0.5237315) q[1];
rz(-pi) q[2];
rz(-0.73170264) q[3];
sx q[3];
rz(-1.9644794) q[3];
sx q[3];
rz(-1.3418152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.61802822) q[2];
sx q[2];
rz(-2.9310493) q[2];
sx q[2];
rz(0.024624126) q[2];
rz(2.5060182) q[3];
sx q[3];
rz(-2.1533951) q[3];
sx q[3];
rz(-2.2583101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(1.4516975) q[0];
sx q[0];
rz(-2.8391333) q[0];
sx q[0];
rz(2.6746993) q[0];
rz(3.0310071) q[1];
sx q[1];
rz(-2.6778335) q[1];
sx q[1];
rz(-1.0708403) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2617874) q[0];
sx q[0];
rz(-2.6222485) q[0];
sx q[0];
rz(-1.0762317) q[0];
rz(0.31857014) q[2];
sx q[2];
rz(-1.397322) q[2];
sx q[2];
rz(-1.8167379) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5731949) q[1];
sx q[1];
rz(-1.0612951) q[1];
sx q[1];
rz(-1.1361109) q[1];
rz(-pi) q[2];
rz(0.91584622) q[3];
sx q[3];
rz(-0.98036843) q[3];
sx q[3];
rz(0.46186033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.11833) q[2];
sx q[2];
rz(-1.7795965) q[2];
sx q[2];
rz(-2.2892717) q[2];
rz(-0.037633745) q[3];
sx q[3];
rz(-1.9084385) q[3];
sx q[3];
rz(-0.5948624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1325876) q[0];
sx q[0];
rz(-1.9629033) q[0];
sx q[0];
rz(-1.8527385) q[0];
rz(-1.5124849) q[1];
sx q[1];
rz(-2.0178724) q[1];
sx q[1];
rz(1.9992794) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90683939) q[0];
sx q[0];
rz(-1.2182353) q[0];
sx q[0];
rz(1.6749188) q[0];
rz(1.8577447) q[2];
sx q[2];
rz(-0.21636886) q[2];
sx q[2];
rz(1.3138101) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3728677) q[1];
sx q[1];
rz(-1.8187675) q[1];
sx q[1];
rz(-2.7889113) q[1];
x q[2];
rz(1.3942424) q[3];
sx q[3];
rz(-1.8774967) q[3];
sx q[3];
rz(-0.97489417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8196572) q[2];
sx q[2];
rz(-1.2677931) q[2];
sx q[2];
rz(0.12040559) q[2];
rz(-1.1558007) q[3];
sx q[3];
rz(-1.6121696) q[3];
sx q[3];
rz(0.22404484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8026546) q[0];
sx q[0];
rz(-1.8010362) q[0];
sx q[0];
rz(-1.4539723) q[0];
rz(1.5441719) q[1];
sx q[1];
rz(-1.6463966) q[1];
sx q[1];
rz(0.45101756) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9266201) q[0];
sx q[0];
rz(-0.30719137) q[0];
sx q[0];
rz(0.19812576) q[0];
rz(-pi) q[1];
rz(2.1233929) q[2];
sx q[2];
rz(-0.42358735) q[2];
sx q[2];
rz(0.46679631) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.34141065) q[1];
sx q[1];
rz(-1.4497611) q[1];
sx q[1];
rz(2.7054663) q[1];
rz(-pi) q[2];
rz(0.0093098442) q[3];
sx q[3];
rz(-2.3311989) q[3];
sx q[3];
rz(1.9805465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1958127) q[2];
sx q[2];
rz(-1.4142298) q[2];
sx q[2];
rz(-1.7808524) q[2];
rz(2.1360548) q[3];
sx q[3];
rz(-1.9660299) q[3];
sx q[3];
rz(1.3895234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0176395) q[0];
sx q[0];
rz(-3.0160976) q[0];
sx q[0];
rz(-0.70607591) q[0];
rz(-1.5438682) q[1];
sx q[1];
rz(-1.4223301) q[1];
sx q[1];
rz(-0.85618883) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5498915) q[0];
sx q[0];
rz(-0.020792637) q[0];
sx q[0];
rz(-1.2077622) q[0];
rz(-pi) q[1];
rz(-1.1107619) q[2];
sx q[2];
rz(-0.80321124) q[2];
sx q[2];
rz(2.3873461) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.013102839) q[1];
sx q[1];
rz(-1.0016514) q[1];
sx q[1];
rz(1.8123958) q[1];
rz(-pi) q[2];
rz(-2.059465) q[3];
sx q[3];
rz(-1.4525982) q[3];
sx q[3];
rz(-1.7080054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84264821) q[2];
sx q[2];
rz(-1.4010669) q[2];
sx q[2];
rz(2.974158) q[2];
rz(-1.5542479) q[3];
sx q[3];
rz(-2.3685679) q[3];
sx q[3];
rz(-1.0003482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-1.7635968) q[0];
sx q[0];
rz(-1.4507699) q[0];
sx q[0];
rz(0.3717306) q[0];
rz(-1.4338214) q[1];
sx q[1];
rz(-1.6038409) q[1];
sx q[1];
rz(-0.30002123) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0473012) q[0];
sx q[0];
rz(-0.30891788) q[0];
sx q[0];
rz(0.39157097) q[0];
rz(-0.062252684) q[2];
sx q[2];
rz(-0.57800284) q[2];
sx q[2];
rz(-2.5684772) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2310861) q[1];
sx q[1];
rz(-0.66201895) q[1];
sx q[1];
rz(0.20259095) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86784243) q[3];
sx q[3];
rz(-1.4757089) q[3];
sx q[3];
rz(-0.33012182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.47664777) q[2];
sx q[2];
rz(-0.46755329) q[2];
sx q[2];
rz(-2.7139968) q[2];
rz(-1.556373) q[3];
sx q[3];
rz(-1.1170324) q[3];
sx q[3];
rz(-2.3868886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72117358) q[0];
sx q[0];
rz(-2.6009646) q[0];
sx q[0];
rz(-2.6721201) q[0];
rz(0.92533127) q[1];
sx q[1];
rz(-1.0877437) q[1];
sx q[1];
rz(-3.1184149) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1184922) q[0];
sx q[0];
rz(-1.7351302) q[0];
sx q[0];
rz(-2.611155) q[0];
rz(-pi) q[1];
x q[1];
rz(2.542664) q[2];
sx q[2];
rz(-2.1492276) q[2];
sx q[2];
rz(-0.79858649) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9564445) q[1];
sx q[1];
rz(-1.7352163) q[1];
sx q[1];
rz(1.7471353) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3599102) q[3];
sx q[3];
rz(-0.74112219) q[3];
sx q[3];
rz(-1.9884584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8699441) q[2];
sx q[2];
rz(-1.9359438) q[2];
sx q[2];
rz(2.6507157) q[2];
rz(2.0513746) q[3];
sx q[3];
rz(-0.96929437) q[3];
sx q[3];
rz(2.8276665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28521095) q[0];
sx q[0];
rz(-0.87775341) q[0];
sx q[0];
rz(-2.0650771) q[0];
rz(2.8648227) q[1];
sx q[1];
rz(-0.77817398) q[1];
sx q[1];
rz(-1.4336817) q[1];
rz(2.4527373) q[2];
sx q[2];
rz(-1.3313455) q[2];
sx q[2];
rz(-0.81305885) q[2];
rz(2.6071045) q[3];
sx q[3];
rz(-1.5535003) q[3];
sx q[3];
rz(1.5737892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
