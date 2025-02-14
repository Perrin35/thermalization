OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.0285398736596107) q[0];
sx q[0];
rz(3.95885035594041) q[0];
sx q[0];
rz(9.13986528515025) q[0];
rz(-2.32183957099915) q[1];
sx q[1];
rz(4.02639624674852) q[1];
sx q[1];
rz(10.5609924554746) q[1];
cx q[1],q[0];
rz(2.32165169715881) q[0];
sx q[0];
rz(3.76430210669572) q[0];
sx q[0];
rz(9.55838150381252) q[0];
rz(-0.89102441072464) q[2];
sx q[2];
rz(5.17007389863069) q[2];
sx q[2];
rz(9.59856869875594) q[2];
cx q[2],q[1];
rz(-0.939282476902008) q[1];
sx q[1];
rz(4.28851810296113) q[1];
sx q[1];
rz(12.2550370454709) q[1];
rz(2.01308631896973) q[3];
sx q[3];
rz(1.77684036095674) q[3];
sx q[3];
rz(8.30863437651798) q[3];
cx q[3],q[2];
rz(3.53743672370911) q[2];
sx q[2];
rz(3.6276006420427) q[2];
sx q[2];
rz(9.13253251313373) q[2];
rz(0.621198356151581) q[3];
sx q[3];
rz(5.20815196831758) q[3];
sx q[3];
rz(9.20602422057792) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.562769651412964) q[0];
sx q[0];
rz(4.25883797009523) q[0];
sx q[0];
rz(9.33836849629089) q[0];
rz(2.74827647209167) q[1];
sx q[1];
rz(4.82912054856355) q[1];
sx q[1];
rz(7.94099185465976) q[1];
cx q[1],q[0];
rz(-0.959796071052551) q[0];
sx q[0];
rz(4.04671839078004) q[0];
sx q[0];
rz(10.6963323116223) q[0];
rz(1.70898354053497) q[2];
sx q[2];
rz(1.88052466710145) q[2];
sx q[2];
rz(10.2781765818517) q[2];
cx q[2],q[1];
rz(0.763689637184143) q[1];
sx q[1];
rz(4.06297067006166) q[1];
sx q[1];
rz(8.04690084456607) q[1];
rz(1.17174053192139) q[3];
sx q[3];
rz(3.74982938368852) q[3];
sx q[3];
rz(7.46892831324741) q[3];
cx q[3],q[2];
rz(0.945744037628174) q[2];
sx q[2];
rz(1.59867790539796) q[2];
sx q[2];
rz(9.42372383900593) q[2];
rz(-0.706831812858582) q[3];
sx q[3];
rz(4.03604188759858) q[3];
sx q[3];
rz(9.70235321520969) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.238667756319046) q[0];
sx q[0];
rz(4.06961235602433) q[0];
sx q[0];
rz(9.903182065479) q[0];
rz(0.844864189624786) q[1];
sx q[1];
rz(4.21790769894654) q[1];
sx q[1];
rz(11.6121430158536) q[1];
cx q[1],q[0];
rz(0.0603928752243519) q[0];
sx q[0];
rz(4.25973740418489) q[0];
sx q[0];
rz(9.84036344885036) q[0];
rz(-0.665047645568848) q[2];
sx q[2];
rz(3.99500647385652) q[2];
sx q[2];
rz(13.1805560350339) q[2];
cx q[2],q[1];
rz(0.409189343452454) q[1];
sx q[1];
rz(5.4848748763376) q[1];
sx q[1];
rz(8.91742250918552) q[1];
rz(-0.266216367483139) q[3];
sx q[3];
rz(2.20532182057435) q[3];
sx q[3];
rz(11.8135323285977) q[3];
cx q[3],q[2];
rz(0.937629520893097) q[2];
sx q[2];
rz(3.48176065285737) q[2];
sx q[2];
rz(9.35350922345325) q[2];
rz(1.01026201248169) q[3];
sx q[3];
rz(4.30260959466035) q[3];
sx q[3];
rz(9.74578935503169) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.456927448511124) q[0];
sx q[0];
rz(4.43428829510743) q[0];
sx q[0];
rz(9.06126765011951) q[0];
rz(-0.892063081264496) q[1];
sx q[1];
rz(1.03282466729219) q[1];
sx q[1];
rz(11.3126417159955) q[1];
cx q[1],q[0];
rz(2.34334683418274) q[0];
sx q[0];
rz(4.32211104233796) q[0];
sx q[0];
rz(9.91637474893733) q[0];
rz(-0.364186525344849) q[2];
sx q[2];
rz(4.61795094807679) q[2];
sx q[2];
rz(10.5254252910535) q[2];
cx q[2],q[1];
rz(2.23347854614258) q[1];
sx q[1];
rz(4.75389054616029) q[1];
sx q[1];
rz(8.10419652461215) q[1];
rz(-0.307800948619843) q[3];
sx q[3];
rz(4.33306899865205) q[3];
sx q[3];
rz(10.0533224701802) q[3];
cx q[3],q[2];
rz(2.75186467170715) q[2];
sx q[2];
rz(3.91565558512742) q[2];
sx q[2];
rz(7.89620540141269) q[2];
rz(-0.450168967247009) q[3];
sx q[3];
rz(4.49377122719819) q[3];
sx q[3];
rz(11.2672909259717) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.124409735202789) q[0];
sx q[0];
rz(2.9876508881622) q[0];
sx q[0];
rz(10.2756567358892) q[0];
rz(1.7393913269043) q[1];
sx q[1];
rz(4.69208541710908) q[1];
sx q[1];
rz(9.57719195484325) q[1];
cx q[1],q[0];
rz(0.352033466100693) q[0];
sx q[0];
rz(4.71855083306367) q[0];
sx q[0];
rz(9.18546337484523) q[0];
rz(0.409884721040726) q[2];
sx q[2];
rz(0.642337950068065) q[2];
sx q[2];
rz(8.32370934485599) q[2];
cx q[2],q[1];
rz(1.72625315189362) q[1];
sx q[1];
rz(4.99498799641664) q[1];
sx q[1];
rz(7.68139991759464) q[1];
rz(0.368666023015976) q[3];
sx q[3];
rz(1.94455018837983) q[3];
sx q[3];
rz(9.71410582064792) q[3];
cx q[3],q[2];
rz(-0.572058141231537) q[2];
sx q[2];
rz(4.06590149004991) q[2];
sx q[2];
rz(9.80064702629253) q[2];
rz(-0.268402606248856) q[3];
sx q[3];
rz(2.12080505688722) q[3];
sx q[3];
rz(11.5275809526364) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.03804504871368) q[0];
sx q[0];
rz(4.05989113648469) q[0];
sx q[0];
rz(9.29902132450744) q[0];
rz(-0.0452995635569096) q[1];
sx q[1];
rz(2.24429753621156) q[1];
sx q[1];
rz(9.92995653151675) q[1];
cx q[1],q[0];
rz(1.40464580059052) q[0];
sx q[0];
rz(2.87897992332513) q[0];
sx q[0];
rz(10.7500247716825) q[0];
rz(-0.491105705499649) q[2];
sx q[2];
rz(4.61507728894288) q[2];
sx q[2];
rz(9.09137657879993) q[2];
cx q[2],q[1];
rz(2.73867011070251) q[1];
sx q[1];
rz(3.94761756260926) q[1];
sx q[1];
rz(9.2218801587741) q[1];
rz(0.466641485691071) q[3];
sx q[3];
rz(1.87100509007508) q[3];
sx q[3];
rz(10.6376125574033) q[3];
cx q[3],q[2];
rz(0.818564116954803) q[2];
sx q[2];
rz(3.47267124255235) q[2];
sx q[2];
rz(6.03007385729953) q[2];
rz(0.744629144668579) q[3];
sx q[3];
rz(4.69534030755097) q[3];
sx q[3];
rz(9.15517727135822) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.44148662686348) q[0];
sx q[0];
rz(4.42811837990815) q[0];
sx q[0];
rz(9.38125335275337) q[0];
rz(1.17319095134735) q[1];
sx q[1];
rz(4.43465343316133) q[1];
sx q[1];
rz(10.2264337301175) q[1];
cx q[1],q[0];
rz(-0.350796818733215) q[0];
sx q[0];
rz(2.96140946646268) q[0];
sx q[0];
rz(10.4136562704961) q[0];
rz(0.84287965297699) q[2];
sx q[2];
rz(4.05090245802934) q[2];
sx q[2];
rz(10.3830295562665) q[2];
cx q[2],q[1];
rz(2.31396269798279) q[1];
sx q[1];
rz(2.45557424624498) q[1];
sx q[1];
rz(9.17825170456573) q[1];
rz(-0.895097374916077) q[3];
sx q[3];
rz(5.90258446534211) q[3];
sx q[3];
rz(12.875698542587) q[3];
cx q[3],q[2];
rz(-0.390485525131226) q[2];
sx q[2];
rz(4.23488000233705) q[2];
sx q[2];
rz(8.50020316838428) q[2];
rz(3.06040024757385) q[3];
sx q[3];
rz(4.8849106152826) q[3];
sx q[3];
rz(10.1528485774915) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.282959491014481) q[0];
sx q[0];
rz(4.38721695740754) q[0];
sx q[0];
rz(10.5002511501233) q[0];
rz(1.76434993743896) q[1];
sx q[1];
rz(4.9110690673166) q[1];
sx q[1];
rz(6.93837258814975) q[1];
cx q[1],q[0];
rz(1.0690883398056) q[0];
sx q[0];
rz(3.07208147098357) q[0];
sx q[0];
rz(10.3149651646535) q[0];
rz(-0.783020079135895) q[2];
sx q[2];
rz(4.50256672699983) q[2];
sx q[2];
rz(10.3736660242002) q[2];
cx q[2],q[1];
rz(0.97444087266922) q[1];
sx q[1];
rz(5.40746107895906) q[1];
sx q[1];
rz(9.33625734447643) q[1];
rz(-0.3386170566082) q[3];
sx q[3];
rz(4.68461814721162) q[3];
sx q[3];
rz(9.59259142576858) q[3];
cx q[3],q[2];
rz(0.367715567350388) q[2];
sx q[2];
rz(3.53312328656251) q[2];
sx q[2];
rz(11.4098402023236) q[2];
rz(1.14415061473846) q[3];
sx q[3];
rz(2.51132276852662) q[3];
sx q[3];
rz(10.7575627326886) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.820412039756775) q[0];
sx q[0];
rz(4.06242874463136) q[0];
sx q[0];
rz(10.6004879236142) q[0];
rz(-1.80204486846924) q[1];
sx q[1];
rz(4.14096984465654) q[1];
sx q[1];
rz(9.83227804898425) q[1];
cx q[1],q[0];
rz(0.202431991696358) q[0];
sx q[0];
rz(3.64939674933488) q[0];
sx q[0];
rz(11.6094197988431) q[0];
rz(-1.0888534784317) q[2];
sx q[2];
rz(2.60662684042985) q[2];
sx q[2];
rz(9.26510140895053) q[2];
cx q[2],q[1];
rz(3.08278512954712) q[1];
sx q[1];
rz(4.73283007939393) q[1];
sx q[1];
rz(10.0123101830403) q[1];
rz(1.84256565570831) q[3];
sx q[3];
rz(4.64376202424104) q[3];
sx q[3];
rz(9.12525517343684) q[3];
cx q[3],q[2];
rz(-1.0743864774704) q[2];
sx q[2];
rz(4.3001379092508) q[2];
sx q[2];
rz(10.4502173423688) q[2];
rz(0.497536718845367) q[3];
sx q[3];
rz(5.3313366492563) q[3];
sx q[3];
rz(9.83124557732745) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.247917264699936) q[0];
sx q[0];
rz(4.29945650895173) q[0];
sx q[0];
rz(11.1809087753217) q[0];
rz(1.04761838912964) q[1];
sx q[1];
rz(4.48649242718751) q[1];
sx q[1];
rz(11.5935322999875) q[1];
cx q[1],q[0];
rz(-0.0948008969426155) q[0];
sx q[0];
rz(3.69229474862153) q[0];
sx q[0];
rz(9.58951922356292) q[0];
rz(-0.734662055969238) q[2];
sx q[2];
rz(3.81807866890962) q[2];
sx q[2];
rz(11.5331926107328) q[2];
cx q[2],q[1];
rz(-0.232158809900284) q[1];
sx q[1];
rz(4.3591834624582) q[1];
sx q[1];
rz(9.57640058397456) q[1];
rz(0.445799469947815) q[3];
sx q[3];
rz(3.8908630331331) q[3];
sx q[3];
rz(9.84749711155101) q[3];
cx q[3],q[2];
rz(-0.732937812805176) q[2];
sx q[2];
rz(3.4698083122545) q[2];
sx q[2];
rz(11.4643842935483) q[2];
rz(1.71056258678436) q[3];
sx q[3];
rz(2.31650629838044) q[3];
sx q[3];
rz(11.3283099889676) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.524845242500305) q[0];
sx q[0];
rz(2.60488137801225) q[0];
sx q[0];
rz(9.21149552463695) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-0.136978030204773) q[1];
sx q[1];
rz(5.15567007859284) q[1];
sx q[1];
rz(10.5272054433744) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.467360764741898) q[2];
sx q[2];
rz(2.0524492581659) q[2];
sx q[2];
rz(9.72134256958171) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.486603885889053) q[3];
sx q[3];
rz(3.48391968210275) q[3];
sx q[3];
rz(8.55506465434238) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
