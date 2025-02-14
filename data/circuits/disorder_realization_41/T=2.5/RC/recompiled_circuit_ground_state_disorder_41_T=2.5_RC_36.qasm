OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1908258) q[0];
sx q[0];
rz(-1.289225) q[0];
sx q[0];
rz(3.0486795) q[0];
rz(-0.095280401) q[1];
sx q[1];
rz(2.4089101) q[1];
sx q[1];
rz(10.664193) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5372972) q[0];
sx q[0];
rz(-1.3193466) q[0];
sx q[0];
rz(2.9049113) q[0];
rz(1.3746512) q[2];
sx q[2];
rz(-1.3117547) q[2];
sx q[2];
rz(-1.0153304) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4463006) q[1];
sx q[1];
rz(-1.9009034) q[1];
sx q[1];
rz(1.6020726) q[1];
rz(-0.69783437) q[3];
sx q[3];
rz(-2.4437332) q[3];
sx q[3];
rz(0.35240155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6713082) q[2];
sx q[2];
rz(-1.4802063) q[2];
sx q[2];
rz(1.3489464) q[2];
rz(0.74364439) q[3];
sx q[3];
rz(-0.27739224) q[3];
sx q[3];
rz(2.9644137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0002366) q[0];
sx q[0];
rz(-1.5970705) q[0];
sx q[0];
rz(2.8479688) q[0];
rz(2.1108421) q[1];
sx q[1];
rz(-1.7799957) q[1];
sx q[1];
rz(0.34293276) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1200877) q[0];
sx q[0];
rz(-1.3653127) q[0];
sx q[0];
rz(1.0387102) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2478825) q[2];
sx q[2];
rz(-2.7512449) q[2];
sx q[2];
rz(-2.5935612) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1520815) q[1];
sx q[1];
rz(-0.73484269) q[1];
sx q[1];
rz(-0.32232018) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9092122) q[3];
sx q[3];
rz(-2.4840601) q[3];
sx q[3];
rz(-0.38207182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.12884101) q[2];
sx q[2];
rz(-1.3920471) q[2];
sx q[2];
rz(-2.3823605) q[2];
rz(-3.1208842) q[3];
sx q[3];
rz(-1.8755251) q[3];
sx q[3];
rz(-2.7032963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.52221209) q[0];
sx q[0];
rz(-0.601957) q[0];
sx q[0];
rz(-3.1371327) q[0];
rz(-2.4941173) q[1];
sx q[1];
rz(-2.5471893) q[1];
sx q[1];
rz(-0.78537816) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68092184) q[0];
sx q[0];
rz(-0.17853949) q[0];
sx q[0];
rz(-2.6723249) q[0];
rz(-pi) q[1];
rz(-1.0713646) q[2];
sx q[2];
rz(-2.7874261) q[2];
sx q[2];
rz(-3.0234697) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8632311) q[1];
sx q[1];
rz(-1.0309217) q[1];
sx q[1];
rz(-3.0616772) q[1];
rz(-pi) q[2];
rz(-2.249751) q[3];
sx q[3];
rz(-1.0247318) q[3];
sx q[3];
rz(-2.8308656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.37041) q[2];
sx q[2];
rz(-0.9202756) q[2];
sx q[2];
rz(1.3578337) q[2];
rz(-0.34058288) q[3];
sx q[3];
rz(-2.1850696) q[3];
sx q[3];
rz(0.79184872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1444645) q[0];
sx q[0];
rz(-2.0576532) q[0];
sx q[0];
rz(2.2564364) q[0];
rz(2.3947233) q[1];
sx q[1];
rz(-2.5179458) q[1];
sx q[1];
rz(1.0964099) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88056662) q[0];
sx q[0];
rz(-1.2231022) q[0];
sx q[0];
rz(-2.8096922) q[0];
rz(-pi) q[1];
x q[1];
rz(1.095216) q[2];
sx q[2];
rz(-3.0747483) q[2];
sx q[2];
rz(-0.16193709) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9982352) q[1];
sx q[1];
rz(-2.0298185) q[1];
sx q[1];
rz(-1.281339) q[1];
rz(-0.73170264) q[3];
sx q[3];
rz(-1.9644794) q[3];
sx q[3];
rz(-1.3418152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61802822) q[2];
sx q[2];
rz(-0.21054331) q[2];
sx q[2];
rz(-3.1169685) q[2];
rz(2.5060182) q[3];
sx q[3];
rz(-2.1533951) q[3];
sx q[3];
rz(-2.2583101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6898952) q[0];
sx q[0];
rz(-2.8391333) q[0];
sx q[0];
rz(2.6746993) q[0];
rz(-0.11058552) q[1];
sx q[1];
rz(-0.46375912) q[1];
sx q[1];
rz(-2.0707524) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.817628) q[0];
sx q[0];
rz(-2.0228798) q[0];
sx q[0];
rz(-2.8766207) q[0];
rz(-pi) q[1];
rz(2.8230225) q[2];
sx q[2];
rz(-1.397322) q[2];
sx q[2];
rz(1.8167379) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9164898) q[1];
sx q[1];
rz(-1.1943294) q[1];
sx q[1];
rz(2.589499) q[1];
rz(2.4038845) q[3];
sx q[3];
rz(-2.2900351) q[3];
sx q[3];
rz(-1.7361452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.11833) q[2];
sx q[2];
rz(-1.3619962) q[2];
sx q[2];
rz(-2.2892717) q[2];
rz(3.1039589) q[3];
sx q[3];
rz(-1.9084385) q[3];
sx q[3];
rz(2.5467303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0090050176) q[0];
sx q[0];
rz(-1.1786893) q[0];
sx q[0];
rz(-1.8527385) q[0];
rz(-1.6291078) q[1];
sx q[1];
rz(-1.1237203) q[1];
sx q[1];
rz(1.9992794) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61296755) q[0];
sx q[0];
rz(-0.36699793) q[0];
sx q[0];
rz(-0.27530833) q[0];
rz(-pi) q[1];
rz(1.283848) q[2];
sx q[2];
rz(-0.21636886) q[2];
sx q[2];
rz(-1.3138101) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9276322) q[1];
sx q[1];
rz(-0.42812706) q[1];
sx q[1];
rz(0.63251782) q[1];
x q[2];
rz(-1.3942424) q[3];
sx q[3];
rz(-1.8774967) q[3];
sx q[3];
rz(0.97489417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8196572) q[2];
sx q[2];
rz(-1.2677931) q[2];
sx q[2];
rz(-3.0211871) q[2];
rz(1.1558007) q[3];
sx q[3];
rz(-1.6121696) q[3];
sx q[3];
rz(2.9175478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33893809) q[0];
sx q[0];
rz(-1.8010362) q[0];
sx q[0];
rz(-1.6876203) q[0];
rz(-1.5441719) q[1];
sx q[1];
rz(-1.495196) q[1];
sx q[1];
rz(0.45101756) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16674834) q[0];
sx q[0];
rz(-1.6303501) q[0];
sx q[0];
rz(-2.8400499) q[0];
rz(-pi) q[1];
rz(-1.2043549) q[2];
sx q[2];
rz(-1.7882573) q[2];
sx q[2];
rz(1.6161473) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9684194) q[1];
sx q[1];
rz(-1.1380769) q[1];
sx q[1];
rz(-1.7041901) q[1];
x q[2];
rz(1.5805832) q[3];
sx q[3];
rz(-0.76044816) q[3];
sx q[3];
rz(1.994054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94577998) q[2];
sx q[2];
rz(-1.7273629) q[2];
sx q[2];
rz(1.7808524) q[2];
rz(2.1360548) q[3];
sx q[3];
rz(-1.9660299) q[3];
sx q[3];
rz(-1.7520693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0176395) q[0];
sx q[0];
rz(-3.0160976) q[0];
sx q[0];
rz(-2.4355167) q[0];
rz(1.5977244) q[1];
sx q[1];
rz(-1.4223301) q[1];
sx q[1];
rz(-0.85618883) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34205758) q[0];
sx q[0];
rz(-1.5634131) q[0];
sx q[0];
rz(1.5513585) q[0];
x q[1];
rz(1.1107619) q[2];
sx q[2];
rz(-0.80321124) q[2];
sx q[2];
rz(-2.3873461) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1284898) q[1];
sx q[1];
rz(-1.0016514) q[1];
sx q[1];
rz(-1.3291969) q[1];
rz(-pi) q[2];
rz(1.3230349) q[3];
sx q[3];
rz(-0.50163402) q[3];
sx q[3];
rz(-0.35546965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.84264821) q[2];
sx q[2];
rz(-1.7405258) q[2];
sx q[2];
rz(-0.16743463) q[2];
rz(1.5873448) q[3];
sx q[3];
rz(-2.3685679) q[3];
sx q[3];
rz(-1.0003482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3779959) q[0];
sx q[0];
rz(-1.4507699) q[0];
sx q[0];
rz(-0.3717306) q[0];
rz(1.7077712) q[1];
sx q[1];
rz(-1.5377518) q[1];
sx q[1];
rz(0.30002123) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094291501) q[0];
sx q[0];
rz(-0.30891788) q[0];
sx q[0];
rz(-2.7500217) q[0];
rz(-pi) q[1];
rz(0.062252684) q[2];
sx q[2];
rz(-2.5635898) q[2];
sx q[2];
rz(0.5731155) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65574291) q[1];
sx q[1];
rz(-2.2169737) q[1];
sx q[1];
rz(-1.7263401) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0172273) q[3];
sx q[3];
rz(-0.87166407) q[3];
sx q[3];
rz(-1.3209526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.47664777) q[2];
sx q[2];
rz(-2.6740394) q[2];
sx q[2];
rz(-2.7139968) q[2];
rz(1.556373) q[3];
sx q[3];
rz(-1.1170324) q[3];
sx q[3];
rz(-0.75470406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72117358) q[0];
sx q[0];
rz(-2.6009646) q[0];
sx q[0];
rz(0.46947259) q[0];
rz(-2.2162614) q[1];
sx q[1];
rz(-2.053849) q[1];
sx q[1];
rz(-0.023177711) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0231004) q[0];
sx q[0];
rz(-1.4064624) q[0];
sx q[0];
rz(-2.611155) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2830354) q[2];
sx q[2];
rz(-2.3344667) q[2];
sx q[2];
rz(-1.6940534) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7851023) q[1];
sx q[1];
rz(-1.7447326) q[1];
sx q[1];
rz(-2.9746303) q[1];
rz(0.18927197) q[3];
sx q[3];
rz(-2.2918275) q[3];
sx q[3];
rz(1.7060351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2716486) q[2];
sx q[2];
rz(-1.2056489) q[2];
sx q[2];
rz(-2.6507157) q[2];
rz(1.0902181) q[3];
sx q[3];
rz(-2.1722983) q[3];
sx q[3];
rz(2.8276665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8563817) q[0];
sx q[0];
rz(-0.87775341) q[0];
sx q[0];
rz(-2.0650771) q[0];
rz(2.8648227) q[1];
sx q[1];
rz(-0.77817398) q[1];
sx q[1];
rz(-1.4336817) q[1];
rz(-1.8770915) q[2];
sx q[2];
rz(-2.2363792) q[2];
sx q[2];
rz(-2.1909942) q[2];
rz(1.5908949) q[3];
sx q[3];
rz(-1.0363967) q[3];
sx q[3];
rz(-3.1283621) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
