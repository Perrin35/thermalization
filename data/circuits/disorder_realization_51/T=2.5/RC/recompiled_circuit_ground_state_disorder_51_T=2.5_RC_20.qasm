OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.92276031) q[0];
sx q[0];
rz(-3.0781167) q[0];
sx q[0];
rz(1.3702962) q[0];
rz(-2.3308347) q[1];
sx q[1];
rz(-1.6914565) q[1];
sx q[1];
rz(2.9210747) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.044174319) q[0];
sx q[0];
rz(-2.3050791) q[0];
sx q[0];
rz(-0.091139779) q[0];
x q[1];
rz(-1.1132672) q[2];
sx q[2];
rz(-1.2731247) q[2];
sx q[2];
rz(0.13099094) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0956362) q[1];
sx q[1];
rz(-1.7929309) q[1];
sx q[1];
rz(1.5820811) q[1];
rz(-pi) q[2];
rz(2.9412782) q[3];
sx q[3];
rz(-1.2216) q[3];
sx q[3];
rz(-1.136245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.084879547) q[2];
sx q[2];
rz(-0.006338174) q[2];
sx q[2];
rz(-2.32178) q[2];
rz(0.020717185) q[3];
sx q[3];
rz(-0.31532225) q[3];
sx q[3];
rz(1.0516385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3725975) q[0];
sx q[0];
rz(-2.9561434) q[0];
sx q[0];
rz(-0.73082596) q[0];
rz(-1.7164879) q[1];
sx q[1];
rz(-2.4162879) q[1];
sx q[1];
rz(2.6740429) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088652164) q[0];
sx q[0];
rz(-0.10693094) q[0];
sx q[0];
rz(-1.8631975) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7317762) q[2];
sx q[2];
rz(-1.6011607) q[2];
sx q[2];
rz(1.8259468) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.0056356965) q[1];
sx q[1];
rz(-1.4219445) q[1];
sx q[1];
rz(0.71133228) q[1];
rz(-pi) q[2];
rz(-2.8465956) q[3];
sx q[3];
rz(-0.85912356) q[3];
sx q[3];
rz(2.0509348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9759489) q[2];
sx q[2];
rz(-3.1303945) q[2];
sx q[2];
rz(2.8051918) q[2];
rz(-0.30763704) q[3];
sx q[3];
rz(-0.010713723) q[3];
sx q[3];
rz(0.78905869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12863185) q[0];
sx q[0];
rz(-0.19879453) q[0];
sx q[0];
rz(3.0096753) q[0];
rz(-1.7031274) q[1];
sx q[1];
rz(-1.4142282) q[1];
sx q[1];
rz(1.6579113) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8298873) q[0];
sx q[0];
rz(-0.97011203) q[0];
sx q[0];
rz(1.2313103) q[0];
rz(1.5533812) q[2];
sx q[2];
rz(-1.546374) q[2];
sx q[2];
rz(-1.2877854) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6710838) q[1];
sx q[1];
rz(-1.6226864) q[1];
sx q[1];
rz(1.4668277) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4999465) q[3];
sx q[3];
rz(-1.6565986) q[3];
sx q[3];
rz(1.5847928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.57054532) q[2];
sx q[2];
rz(-1.3344301) q[2];
sx q[2];
rz(-1.6837616) q[2];
rz(0.7782065) q[3];
sx q[3];
rz(-3.1359378) q[3];
sx q[3];
rz(0.2125423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.0197765) q[0];
sx q[0];
rz(-1.7541405) q[0];
sx q[0];
rz(-0.29086599) q[0];
rz(1.5485171) q[1];
sx q[1];
rz(-0.80268186) q[1];
sx q[1];
rz(-3.1150418) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8940116) q[0];
sx q[0];
rz(-2.6373626) q[0];
sx q[0];
rz(-0.79968851) q[0];
x q[1];
rz(-0.45088972) q[2];
sx q[2];
rz(-1.5227382) q[2];
sx q[2];
rz(-2.5250736) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2902879) q[1];
sx q[1];
rz(-1.0920224) q[1];
sx q[1];
rz(-2.4393875) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7351355) q[3];
sx q[3];
rz(-3.1222635) q[3];
sx q[3];
rz(0.011474495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.87963858) q[2];
sx q[2];
rz(-0.018435437) q[2];
sx q[2];
rz(-2.7035942) q[2];
rz(-1.1079463) q[3];
sx q[3];
rz(-3.1236533) q[3];
sx q[3];
rz(-1.3560791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9172025) q[0];
sx q[0];
rz(-2.6438535) q[0];
sx q[0];
rz(-2.9955731) q[0];
rz(-3.0085425) q[1];
sx q[1];
rz(-1.0597884) q[1];
sx q[1];
rz(2.6859247) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9216292) q[0];
sx q[0];
rz(-1.5507924) q[0];
sx q[0];
rz(3.0476493) q[0];
rz(-0.44134042) q[2];
sx q[2];
rz(-0.63950634) q[2];
sx q[2];
rz(-0.17989448) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7130374) q[1];
sx q[1];
rz(-0.26838612) q[1];
sx q[1];
rz(2.1191855) q[1];
rz(2.9140329) q[3];
sx q[3];
rz(-1.5557845) q[3];
sx q[3];
rz(0.14666569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2666152) q[2];
sx q[2];
rz(-3.1217323) q[2];
sx q[2];
rz(-1.88545) q[2];
rz(0.52136326) q[3];
sx q[3];
rz(-0.25786906) q[3];
sx q[3];
rz(-0.34463394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94619036) q[0];
sx q[0];
rz(-2.7552216) q[0];
sx q[0];
rz(-2.572701) q[0];
rz(-0.16828123) q[1];
sx q[1];
rz(-1.556309) q[1];
sx q[1];
rz(3.1313484) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55355799) q[0];
sx q[0];
rz(-1.7419683) q[0];
sx q[0];
rz(-1.5215015) q[0];
rz(-pi) q[1];
rz(3.1222759) q[2];
sx q[2];
rz(-1.7512808) q[2];
sx q[2];
rz(1.7648802) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3553487) q[1];
sx q[1];
rz(-3.0867483) q[1];
sx q[1];
rz(2.2516903) q[1];
x q[2];
rz(-1.5404352) q[3];
sx q[3];
rz(-0.93685461) q[3];
sx q[3];
rz(0.46677865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9804618) q[2];
sx q[2];
rz(-0.22393301) q[2];
sx q[2];
rz(1.0978318) q[2];
rz(0.3499507) q[3];
sx q[3];
rz(-1.8643458) q[3];
sx q[3];
rz(2.2866975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.8769281) q[0];
sx q[0];
rz(-2.9783037) q[0];
sx q[0];
rz(-2.7352585) q[0];
rz(1.8304652) q[1];
sx q[1];
rz(-0.51478148) q[1];
sx q[1];
rz(-3.0313671) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12719181) q[0];
sx q[0];
rz(-1.4716291) q[0];
sx q[0];
rz(-1.5983461) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1821457) q[2];
sx q[2];
rz(-3.1350208) q[2];
sx q[2];
rz(-0.98192865) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4496042) q[1];
sx q[1];
rz(-1.1432198) q[1];
sx q[1];
rz(3.1319992) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87894999) q[3];
sx q[3];
rz(-0.96262299) q[3];
sx q[3];
rz(-1.5076306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8608287) q[2];
sx q[2];
rz(-3.132756) q[2];
sx q[2];
rz(-2.3542118) q[2];
rz(-0.27960882) q[3];
sx q[3];
rz(-3.0968554) q[3];
sx q[3];
rz(2.4217822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-3.1396007) q[0];
sx q[0];
rz(-0.75374341) q[0];
sx q[0];
rz(1.7107704) q[0];
rz(0.17499533) q[1];
sx q[1];
rz(-2.4038959) q[1];
sx q[1];
rz(1.4281248) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.437182) q[0];
sx q[0];
rz(-1.5927654) q[0];
sx q[0];
rz(3.1385413) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0753209) q[2];
sx q[2];
rz(-3.1331535) q[2];
sx q[2];
rz(-1.5167459) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0261532) q[1];
sx q[1];
rz(-0.31603957) q[1];
sx q[1];
rz(-2.1137966) q[1];
x q[2];
rz(1.346502) q[3];
sx q[3];
rz(-1.845775) q[3];
sx q[3];
rz(-1.1077653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3565107) q[2];
sx q[2];
rz(-0.76866895) q[2];
sx q[2];
rz(1.53995) q[2];
rz(-2.8421863) q[3];
sx q[3];
rz(-1.4842024) q[3];
sx q[3];
rz(-0.33777344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47141075) q[0];
sx q[0];
rz(-0.010936745) q[0];
sx q[0];
rz(-2.6543044) q[0];
rz(-1.7022853) q[1];
sx q[1];
rz(-0.7936365) q[1];
sx q[1];
rz(-2.7557441) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4131933) q[0];
sx q[0];
rz(-1.5110713) q[0];
sx q[0];
rz(2.7483536) q[0];
rz(-1.9963485) q[2];
sx q[2];
rz(-1.4686079) q[2];
sx q[2];
rz(-2.9152753) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0688975) q[1];
sx q[1];
rz(-2.4494315) q[1];
sx q[1];
rz(0.051874224) q[1];
rz(-pi) q[2];
rz(-3.039592) q[3];
sx q[3];
rz(-2.8047562) q[3];
sx q[3];
rz(-0.61397314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4243329) q[2];
sx q[2];
rz(-0.029742664) q[2];
sx q[2];
rz(-2.7359803) q[2];
rz(-0.82617104) q[3];
sx q[3];
rz(-3.0982389) q[3];
sx q[3];
rz(-0.9894754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9721603) q[0];
sx q[0];
rz(-2.5763474) q[0];
sx q[0];
rz(-1.3954847) q[0];
rz(-1.6204429) q[1];
sx q[1];
rz(-0.65137678) q[1];
sx q[1];
rz(1.3491389) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5842218) q[0];
sx q[0];
rz(-2.2855846) q[0];
sx q[0];
rz(2.8609311) q[0];
x q[1];
rz(-1.8113891) q[2];
sx q[2];
rz(-1.7363915) q[2];
sx q[2];
rz(0.22364932) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.031489109) q[1];
sx q[1];
rz(-1.5851136) q[1];
sx q[1];
rz(-1.5107461) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1114499) q[3];
sx q[3];
rz(-1.5858012) q[3];
sx q[3];
rz(2.4001013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.32538432) q[2];
sx q[2];
rz(-3.1367229) q[2];
sx q[2];
rz(1.4332786) q[2];
rz(1.2895182) q[3];
sx q[3];
rz(-3.0472445) q[3];
sx q[3];
rz(-1.2963699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0636487) q[0];
sx q[0];
rz(-2.1049121) q[0];
sx q[0];
rz(-2.5459469) q[0];
rz(3.0972277) q[1];
sx q[1];
rz(-2.8652419) q[1];
sx q[1];
rz(-1.2038632) q[1];
rz(-1.9738214) q[2];
sx q[2];
rz(-2.6296305) q[2];
sx q[2];
rz(0.59811022) q[2];
rz(1.1607004) q[3];
sx q[3];
rz(-1.6312508) q[3];
sx q[3];
rz(-1.4283258) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
