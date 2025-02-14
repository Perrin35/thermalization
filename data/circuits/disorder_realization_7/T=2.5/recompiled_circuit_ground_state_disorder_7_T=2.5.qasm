OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.112438417971134) q[0];
sx q[0];
rz(4.7726790030771) q[0];
sx q[0];
rz(9.79092380999728) q[0];
rz(2.08682990074158) q[1];
sx q[1];
rz(4.42594054539735) q[1];
sx q[1];
rz(8.68518671988651) q[1];
cx q[1],q[0];
rz(-0.306924313306808) q[0];
sx q[0];
rz(3.73605492909486) q[0];
sx q[0];
rz(10.0905762076299) q[0];
rz(-0.934112727642059) q[2];
sx q[2];
rz(4.21239856083924) q[2];
sx q[2];
rz(10.6186597108762) q[2];
cx q[2],q[1];
rz(-2.22409868240356) q[1];
sx q[1];
rz(1.98278942902619) q[1];
sx q[1];
rz(12.0068406820218) q[1];
rz(0.876189231872559) q[3];
sx q[3];
rz(2.73121535976464) q[3];
sx q[3];
rz(9.25113894640609) q[3];
cx q[3],q[2];
rz(0.710204839706421) q[2];
sx q[2];
rz(4.57916918595368) q[2];
sx q[2];
rz(10.5064964055936) q[2];
rz(1.2181351184845) q[3];
sx q[3];
rz(4.7033248265558) q[3];
sx q[3];
rz(10.8959878444593) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.00649380683899) q[0];
sx q[0];
rz(2.30689373810823) q[0];
sx q[0];
rz(10.2619686484258) q[0];
rz(0.937650144100189) q[1];
sx q[1];
rz(4.52719393570954) q[1];
sx q[1];
rz(9.29991295038863) q[1];
cx q[1],q[0];
rz(2.36280846595764) q[0];
sx q[0];
rz(1.15114012558992) q[0];
sx q[0];
rz(10.4383505344312) q[0];
rz(4.26804113388062) q[2];
sx q[2];
rz(3.76717922289903) q[2];
sx q[2];
rz(4.56688878535434) q[2];
cx q[2],q[1];
rz(-1.30954504013062) q[1];
sx q[1];
rz(4.62859133084352) q[1];
sx q[1];
rz(11.5970759153287) q[1];
rz(1.45785439014435) q[3];
sx q[3];
rz(3.4491502066427) q[3];
sx q[3];
rz(8.1764045715253) q[3];
cx q[3],q[2];
rz(-0.798959255218506) q[2];
sx q[2];
rz(5.21710935433442) q[2];
sx q[2];
rz(9.01808977722331) q[2];
rz(-1.13253557682037) q[3];
sx q[3];
rz(4.66038289864595) q[3];
sx q[3];
rz(10.2922160983007) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.366585284471512) q[0];
sx q[0];
rz(3.72713384230668) q[0];
sx q[0];
rz(10.5051471948545) q[0];
rz(-0.222515463829041) q[1];
sx q[1];
rz(2.37189260323579) q[1];
sx q[1];
rz(9.99535987376376) q[1];
cx q[1],q[0];
rz(1.27043175697327) q[0];
sx q[0];
rz(5.78245416482026) q[0];
sx q[0];
rz(11.8807835340421) q[0];
rz(-2.16574454307556) q[2];
sx q[2];
rz(5.24681940873201) q[2];
sx q[2];
rz(9.79402939080402) q[2];
cx q[2],q[1];
rz(0.966484427452087) q[1];
sx q[1];
rz(5.47982588608796) q[1];
sx q[1];
rz(11.0262114763181) q[1];
rz(-0.181811943650246) q[3];
sx q[3];
rz(5.27985730965669) q[3];
sx q[3];
rz(9.79566363095447) q[3];
cx q[3],q[2];
rz(-0.835016191005707) q[2];
sx q[2];
rz(4.62275782425935) q[2];
sx q[2];
rz(11.331462955467) q[2];
rz(1.07704770565033) q[3];
sx q[3];
rz(4.70051589806611) q[3];
sx q[3];
rz(9.00748679637119) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.00175891048274934) q[0];
sx q[0];
rz(4.91248324711854) q[0];
sx q[0];
rz(11.8860780954282) q[0];
rz(1.74400496482849) q[1];
sx q[1];
rz(2.00605073769624) q[1];
sx q[1];
rz(9.18955457805797) q[1];
cx q[1],q[0];
rz(-0.096073791384697) q[0];
sx q[0];
rz(4.17291286786134) q[0];
sx q[0];
rz(9.84811297654315) q[0];
rz(-0.716032326221466) q[2];
sx q[2];
rz(1.82734421093995) q[2];
sx q[2];
rz(13.5046858549039) q[2];
cx q[2],q[1];
rz(-0.669162333011627) q[1];
sx q[1];
rz(4.14411965210969) q[1];
sx q[1];
rz(9.16394043564006) q[1];
rz(1.19182193279266) q[3];
sx q[3];
rz(6.56977740128572) q[3];
sx q[3];
rz(9.43091358932807) q[3];
cx q[3],q[2];
rz(0.325171649456024) q[2];
sx q[2];
rz(5.20724526246125) q[2];
sx q[2];
rz(8.04421827792331) q[2];
rz(0.766624629497528) q[3];
sx q[3];
rz(3.86150965292985) q[3];
sx q[3];
rz(10.0648032188336) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.04708731174469) q[0];
sx q[0];
rz(2.6386058648401) q[0];
sx q[0];
rz(10.692730998985) q[0];
rz(0.937159299850464) q[1];
sx q[1];
rz(4.07633188565309) q[1];
sx q[1];
rz(9.68235010503932) q[1];
cx q[1],q[0];
rz(0.661438882350922) q[0];
sx q[0];
rz(4.15697279770906) q[0];
sx q[0];
rz(11.1737720727841) q[0];
rz(4.40400314331055) q[2];
sx q[2];
rz(4.13333043654496) q[2];
sx q[2];
rz(5.90570614337131) q[2];
cx q[2],q[1];
rz(-1.39933812618256) q[1];
sx q[1];
rz(5.78708520730073) q[1];
sx q[1];
rz(10.5401302337567) q[1];
rz(0.42970284819603) q[3];
sx q[3];
rz(2.61670157511766) q[3];
sx q[3];
rz(10.4683460950772) q[3];
cx q[3],q[2];
rz(-0.816897213459015) q[2];
sx q[2];
rz(2.2163052876764) q[2];
sx q[2];
rz(10.2683601140897) q[2];
rz(0.905136108398438) q[3];
sx q[3];
rz(4.03750452597673) q[3];
sx q[3];
rz(9.68824333547756) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.05966854095459) q[0];
sx q[0];
rz(2.79842394788797) q[0];
sx q[0];
rz(10.8610197067182) q[0];
rz(-1.24172711372375) q[1];
sx q[1];
rz(1.42715862591798) q[1];
sx q[1];
rz(11.9480755090634) q[1];
cx q[1],q[0];
rz(0.70253187417984) q[0];
sx q[0];
rz(2.09053221543367) q[0];
sx q[0];
rz(8.48829040526553) q[0];
rz(1.73617994785309) q[2];
sx q[2];
rz(4.03674665291841) q[2];
sx q[2];
rz(9.84212905763789) q[2];
cx q[2],q[1];
rz(0.306963175535202) q[1];
sx q[1];
rz(4.41119721730287) q[1];
sx q[1];
rz(8.28338274954959) q[1];
rz(1.07384729385376) q[3];
sx q[3];
rz(2.59223106701905) q[3];
sx q[3];
rz(9.55775949954196) q[3];
cx q[3],q[2];
rz(0.618784070014954) q[2];
sx q[2];
rz(4.38713100750978) q[2];
sx q[2];
rz(10.6194321870725) q[2];
rz(-1.49058020114899) q[3];
sx q[3];
rz(4.5574690421396) q[3];
sx q[3];
rz(8.13344631194278) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.24565827846527) q[0];
sx q[0];
rz(2.77225616772706) q[0];
sx q[0];
rz(9.58822993039295) q[0];
rz(0.555119037628174) q[1];
sx q[1];
rz(4.6149469931894) q[1];
sx q[1];
rz(8.60094264744922) q[1];
cx q[1],q[0];
rz(1.46519827842712) q[0];
sx q[0];
rz(4.85611096222932) q[0];
sx q[0];
rz(11.3149622440259) q[0];
rz(1.83748507499695) q[2];
sx q[2];
rz(5.47754612763459) q[2];
sx q[2];
rz(11.4476890325467) q[2];
cx q[2],q[1];
rz(-0.403856217861176) q[1];
sx q[1];
rz(4.25940934022004) q[1];
sx q[1];
rz(9.62286905049487) q[1];
rz(-0.904558122158051) q[3];
sx q[3];
rz(3.64842072327668) q[3];
sx q[3];
rz(12.6586568117063) q[3];
cx q[3],q[2];
rz(1.05503630638123) q[2];
sx q[2];
rz(4.44724372227723) q[2];
sx q[2];
rz(10.8894642352979) q[2];
rz(3.08754086494446) q[3];
sx q[3];
rz(2.26501044829423) q[3];
sx q[3];
rz(9.00211060642406) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.61580491065979) q[0];
sx q[0];
rz(3.52989649971063) q[0];
sx q[0];
rz(10.8565100192945) q[0];
rz(1.88523852825165) q[1];
sx q[1];
rz(4.59952977498109) q[1];
sx q[1];
rz(11.2938301324765) q[1];
cx q[1],q[0];
rz(-1.15210688114166) q[0];
sx q[0];
rz(2.47403344710404) q[0];
sx q[0];
rz(10.7555489301602) q[0];
rz(-0.283207625150681) q[2];
sx q[2];
rz(4.22134867508943) q[2];
sx q[2];
rz(10.5839516878049) q[2];
cx q[2],q[1];
rz(1.14507293701172) q[1];
sx q[1];
rz(3.61218348343904) q[1];
sx q[1];
rz(8.63604513405963) q[1];
rz(-0.813855350017548) q[3];
sx q[3];
rz(4.45234862168367) q[3];
sx q[3];
rz(9.50323749183818) q[3];
cx q[3],q[2];
rz(0.584974765777588) q[2];
sx q[2];
rz(1.40404001076753) q[2];
sx q[2];
rz(7.00967357157871) q[2];
rz(-0.118974208831787) q[3];
sx q[3];
rz(4.53468218644197) q[3];
sx q[3];
rz(9.94550648926898) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.93100905418396) q[0];
sx q[0];
rz(1.98876610596711) q[0];
sx q[0];
rz(9.7782901585023) q[0];
rz(1.17342364788055) q[1];
sx q[1];
rz(4.42720809777314) q[1];
sx q[1];
rz(8.13982770442172) q[1];
cx q[1],q[0];
rz(1.37473201751709) q[0];
sx q[0];
rz(3.59997284610803) q[0];
sx q[0];
rz(9.10920736788913) q[0];
rz(-0.859734296798706) q[2];
sx q[2];
rz(3.71104327042634) q[2];
sx q[2];
rz(10.5013557434003) q[2];
cx q[2],q[1];
rz(1.35815298557281) q[1];
sx q[1];
rz(3.43647426565225) q[1];
sx q[1];
rz(9.35031843035623) q[1];
rz(1.33416509628296) q[3];
sx q[3];
rz(3.77156481345231) q[3];
sx q[3];
rz(9.55554642378494) q[3];
cx q[3],q[2];
rz(0.533591985702515) q[2];
sx q[2];
rz(3.9777462204271) q[2];
sx q[2];
rz(11.7509214639585) q[2];
rz(0.647255301475525) q[3];
sx q[3];
rz(4.904689462977) q[3];
sx q[3];
rz(10.2063981652181) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.6668928861618) q[0];
sx q[0];
rz(6.01673069794709) q[0];
sx q[0];
rz(11.1347743034284) q[0];
rz(3.55309534072876) q[1];
sx q[1];
rz(5.13463059266145) q[1];
sx q[1];
rz(9.55327548681899) q[1];
cx q[1],q[0];
rz(-2.08138251304626) q[0];
sx q[0];
rz(1.8888705094629) q[0];
sx q[0];
rz(10.312368130676) q[0];
rz(-2.43120121955872) q[2];
sx q[2];
rz(1.67710724671418) q[2];
sx q[2];
rz(12.8436903715055) q[2];
cx q[2],q[1];
rz(0.656622767448425) q[1];
sx q[1];
rz(4.04106405575807) q[1];
sx q[1];
rz(10.0976434707563) q[1];
rz(-0.413025200366974) q[3];
sx q[3];
rz(4.22261956532533) q[3];
sx q[3];
rz(9.86492172478839) q[3];
cx q[3],q[2];
rz(-1.07820415496826) q[2];
sx q[2];
rz(4.95606890519197) q[2];
sx q[2];
rz(8.64588842391177) q[2];
rz(1.28634595870972) q[3];
sx q[3];
rz(4.23639336426789) q[3];
sx q[3];
rz(10.6902279615323) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.97958266735077) q[0];
sx q[0];
rz(2.73382863600785) q[0];
sx q[0];
rz(10.4710396289746) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-1.04142487049103) q[1];
sx q[1];
rz(2.61905417044694) q[1];
sx q[1];
rz(9.59228738247558) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.731695652008057) q[2];
sx q[2];
rz(3.70377013285691) q[2];
sx q[2];
rz(11.7486645936887) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.783587574958801) q[3];
sx q[3];
rz(2.33681485255296) q[3];
sx q[3];
rz(9.48533252476856) q[3];
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
