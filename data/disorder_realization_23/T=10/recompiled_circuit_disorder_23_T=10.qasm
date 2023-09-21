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
rz(1.26205790042877) q[0];
sx q[0];
rz(7.69277778466279) q[0];
sx q[0];
rz(11.1322439670484) q[0];
rz(0.634245097637177) q[1];
sx q[1];
rz(6.8847817500406) q[1];
sx q[1];
rz(9.84319815634891) q[1];
cx q[1],q[0];
rz(2.3474862575531) q[0];
sx q[0];
rz(5.73880902131135) q[0];
sx q[0];
rz(8.8988406419675) q[0];
rz(1.68249177932739) q[2];
sx q[2];
rz(1.70614913304383) q[2];
sx q[2];
rz(6.68776724337741) q[2];
cx q[2],q[1];
rz(0.132222086191177) q[1];
sx q[1];
rz(1.81440606911714) q[1];
sx q[1];
rz(12.0974957704465) q[1];
rz(3.58742499351501) q[3];
sx q[3];
rz(3.62164458830888) q[3];
sx q[3];
rz(8.30406937598392) q[3];
cx q[3],q[2];
rz(1.81995403766632) q[2];
sx q[2];
rz(4.51077643235261) q[2];
sx q[2];
rz(11.7283951997678) q[2];
rz(-0.493012458086014) q[3];
sx q[3];
rz(3.41450706322724) q[3];
sx q[3];
rz(9.5037698879759) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.747196018695831) q[0];
sx q[0];
rz(2.41737392743165) q[0];
sx q[0];
rz(10.7026285886685) q[0];
rz(-0.176780745387077) q[1];
sx q[1];
rz(4.96880200703675) q[1];
sx q[1];
rz(6.71535251139804) q[1];
cx q[1],q[0];
rz(-1.99581873416901) q[0];
sx q[0];
rz(5.6818951686197) q[0];
sx q[0];
rz(8.35651252268955) q[0];
rz(1.17064332962036) q[2];
sx q[2];
rz(1.17409852345521) q[2];
sx q[2];
rz(6.47361896037265) q[2];
cx q[2],q[1];
rz(1.19516158103943) q[1];
sx q[1];
rz(5.8257533629709) q[1];
sx q[1];
rz(12.2990090608518) q[1];
rz(0.565000832080841) q[3];
sx q[3];
rz(4.73338619072969) q[3];
sx q[3];
rz(10.2526082754056) q[3];
cx q[3],q[2];
rz(-0.549234867095947) q[2];
sx q[2];
rz(1.26404145558412) q[2];
sx q[2];
rz(12.0796706437986) q[2];
rz(-1.37820792198181) q[3];
sx q[3];
rz(4.40156522591645) q[3];
sx q[3];
rz(8.89195582865878) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0407226718962193) q[0];
sx q[0];
rz(3.88807985384996) q[0];
sx q[0];
rz(9.84212431906863) q[0];
rz(4.63025712966919) q[1];
sx q[1];
rz(3.68709102471406) q[1];
sx q[1];
rz(6.78957030772373) q[1];
cx q[1],q[0];
rz(-0.157841414213181) q[0];
sx q[0];
rz(6.67960062821443) q[0];
sx q[0];
rz(8.95197788476154) q[0];
rz(2.30531001091003) q[2];
sx q[2];
rz(1.75064149697358) q[2];
sx q[2];
rz(10.1339245796125) q[2];
cx q[2],q[1];
rz(-3.03269290924072) q[1];
sx q[1];
rz(4.27321317990357) q[1];
sx q[1];
rz(9.76591250895664) q[1];
rz(-0.369549959897995) q[3];
sx q[3];
rz(4.54196420510346) q[3];
sx q[3];
rz(8.50951150654956) q[3];
cx q[3],q[2];
rz(-0.699040770530701) q[2];
sx q[2];
rz(3.60295084317262) q[2];
sx q[2];
rz(8.83330527543231) q[2];
rz(0.586020350456238) q[3];
sx q[3];
rz(1.20892170270021) q[3];
sx q[3];
rz(11.135192012779) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.16990029811859) q[0];
sx q[0];
rz(2.51328918536241) q[0];
sx q[0];
rz(10.263961172096) q[0];
rz(-3.16717123985291) q[1];
sx q[1];
rz(0.695681007700511) q[1];
sx q[1];
rz(7.87619409560367) q[1];
cx q[1],q[0];
rz(2.4913113117218) q[0];
sx q[0];
rz(2.68227526743943) q[0];
sx q[0];
rz(9.21828704177543) q[0];
rz(0.726857542991638) q[2];
sx q[2];
rz(5.37043753464753) q[2];
sx q[2];
rz(9.07998535632297) q[2];
cx q[2],q[1];
rz(-0.384944260120392) q[1];
sx q[1];
rz(5.09458366234834) q[1];
sx q[1];
rz(10.2274759769361) q[1];
rz(1.14021301269531) q[3];
sx q[3];
rz(4.97084823449189) q[3];
sx q[3];
rz(11.9900977373044) q[3];
cx q[3],q[2];
rz(2.83623909950256) q[2];
sx q[2];
rz(5.39919415314729) q[2];
sx q[2];
rz(9.3250918969433) q[2];
rz(2.18274068832397) q[3];
sx q[3];
rz(1.8226262648874) q[3];
sx q[3];
rz(8.09980390071079) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.56617671251297) q[0];
sx q[0];
rz(1.75949362118775) q[0];
sx q[0];
rz(9.16883755325481) q[0];
rz(2.68049645423889) q[1];
sx q[1];
rz(5.23950424988801) q[1];
sx q[1];
rz(10.1848392844121) q[1];
cx q[1],q[0];
rz(1.59986591339111) q[0];
sx q[0];
rz(4.59562245209748) q[0];
sx q[0];
rz(9.52852411418363) q[0];
rz(2.72694325447083) q[2];
sx q[2];
rz(4.65412512620027) q[2];
sx q[2];
rz(10.525365805618) q[2];
cx q[2],q[1];
rz(0.45436817407608) q[1];
sx q[1];
rz(4.83905795414979) q[1];
sx q[1];
rz(9.29578063487216) q[1];
rz(-1.93645536899567) q[3];
sx q[3];
rz(5.32620373566682) q[3];
sx q[3];
rz(10.3957429885785) q[3];
cx q[3],q[2];
rz(-0.570506751537323) q[2];
sx q[2];
rz(4.91082290013368) q[2];
sx q[2];
rz(8.75053826569721) q[2];
rz(0.214806020259857) q[3];
sx q[3];
rz(2.68476969202096) q[3];
sx q[3];
rz(9.44212230703934) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.531334161758423) q[0];
sx q[0];
rz(4.81275442441041) q[0];
sx q[0];
rz(11.387256717674) q[0];
rz(0.204826518893242) q[1];
sx q[1];
rz(2.34634593327577) q[1];
sx q[1];
rz(11.4994101285855) q[1];
cx q[1],q[0];
rz(1.24440062046051) q[0];
sx q[0];
rz(3.16684164677794) q[0];
sx q[0];
rz(10.3841570377271) q[0];
rz(1.14588057994843) q[2];
sx q[2];
rz(4.30853870709474) q[2];
sx q[2];
rz(7.34769604205295) q[2];
cx q[2],q[1];
rz(-1.30315911769867) q[1];
sx q[1];
rz(4.66488829453523) q[1];
sx q[1];
rz(10.171931898586) q[1];
rz(0.163192763924599) q[3];
sx q[3];
rz(2.36051514943177) q[3];
sx q[3];
rz(10.6738381147306) q[3];
cx q[3],q[2];
rz(-0.710104644298553) q[2];
sx q[2];
rz(4.4308638890558) q[2];
sx q[2];
rz(9.98472121953174) q[2];
rz(2.41527485847473) q[3];
sx q[3];
rz(2.83282187779481) q[3];
sx q[3];
rz(9.7302755176942) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.85754919052124) q[0];
sx q[0];
rz(3.68815848429734) q[0];
sx q[0];
rz(8.00366792678043) q[0];
rz(0.202060848474503) q[1];
sx q[1];
rz(4.57544902165467) q[1];
sx q[1];
rz(10.2829546689908) q[1];
cx q[1],q[0];
rz(1.21505427360535) q[0];
sx q[0];
rz(3.33715307910974) q[0];
sx q[0];
rz(10.1441263914029) q[0];
rz(0.103262327611446) q[2];
sx q[2];
rz(2.72913602192933) q[2];
sx q[2];
rz(9.02771372198268) q[2];
cx q[2],q[1];
rz(4.89942169189453) q[1];
sx q[1];
rz(4.74135628541047) q[1];
sx q[1];
rz(9.44005513693347) q[1];
rz(0.34595975279808) q[3];
sx q[3];
rz(5.25609532197053) q[3];
sx q[3];
rz(12.916548228256) q[3];
cx q[3],q[2];
rz(0.287856072187424) q[2];
sx q[2];
rz(2.66174721916253) q[2];
sx q[2];
rz(10.7502114534299) q[2];
rz(0.89007967710495) q[3];
sx q[3];
rz(4.28869405587251) q[3];
sx q[3];
rz(8.4362551331441) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.67780745029449) q[0];
sx q[0];
rz(3.71414187748963) q[0];
sx q[0];
rz(9.05006397365733) q[0];
rz(-2.1627140045166) q[1];
sx q[1];
rz(0.681904943781444) q[1];
sx q[1];
rz(11.2168645620267) q[1];
cx q[1],q[0];
rz(-0.230962023139) q[0];
sx q[0];
rz(2.17191425164277) q[0];
sx q[0];
rz(9.89887673257991) q[0];
rz(1.5970686674118) q[2];
sx q[2];
rz(0.718258293467112) q[2];
sx q[2];
rz(12.7657019853513) q[2];
cx q[2],q[1];
rz(3.46019625663757) q[1];
sx q[1];
rz(3.33012569149072) q[1];
sx q[1];
rz(9.1573242008607) q[1];
rz(-1.98614454269409) q[3];
sx q[3];
rz(6.8691715319925) q[3];
sx q[3];
rz(10.677022433273) q[3];
cx q[3],q[2];
rz(2.49020409584045) q[2];
sx q[2];
rz(5.53043285210664) q[2];
sx q[2];
rz(12.0976748227994) q[2];
rz(1.19413411617279) q[3];
sx q[3];
rz(4.4457302411371) q[3];
sx q[3];
rz(8.06118545531436) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.630121827125549) q[0];
sx q[0];
rz(4.04793581564958) q[0];
sx q[0];
rz(11.304467535011) q[0];
rz(0.175038605928421) q[1];
sx q[1];
rz(5.14134541352326) q[1];
sx q[1];
rz(10.9623465299527) q[1];
cx q[1],q[0];
rz(-0.810485482215881) q[0];
sx q[0];
rz(5.01663485367829) q[0];
sx q[0];
rz(10.7637609004895) q[0];
rz(1.34013915061951) q[2];
sx q[2];
rz(2.71291774709756) q[2];
sx q[2];
rz(5.9603583574216) q[2];
cx q[2],q[1];
rz(1.69198417663574) q[1];
sx q[1];
rz(2.27126279671723) q[1];
sx q[1];
rz(9.57335708140537) q[1];
rz(0.650534331798553) q[3];
sx q[3];
rz(2.58396211464936) q[3];
sx q[3];
rz(9.87093958853885) q[3];
cx q[3],q[2];
rz(0.0399154461920261) q[2];
sx q[2];
rz(4.09414443572099) q[2];
sx q[2];
rz(10.1861227512281) q[2];
rz(0.900416553020477) q[3];
sx q[3];
rz(5.68368545373017) q[3];
sx q[3];
rz(9.47380193918153) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.339235663414) q[0];
sx q[0];
rz(3.60986331303651) q[0];
sx q[0];
rz(9.20787217318221) q[0];
rz(0.631981730461121) q[1];
sx q[1];
rz(4.63302996953065) q[1];
sx q[1];
rz(8.47004715203449) q[1];
cx q[1],q[0];
rz(-1.27495980262756) q[0];
sx q[0];
rz(5.0325517972284) q[0];
sx q[0];
rz(10.2352072358052) q[0];
rz(-0.940073609352112) q[2];
sx q[2];
rz(5.63436833222444) q[2];
sx q[2];
rz(11.4479031324308) q[2];
cx q[2],q[1];
rz(0.732505142688751) q[1];
sx q[1];
rz(5.1675023158365) q[1];
sx q[1];
rz(9.50959729253455) q[1];
rz(-1.12761890888214) q[3];
sx q[3];
rz(4.12583235104615) q[3];
sx q[3];
rz(9.29541111587688) q[3];
cx q[3],q[2];
rz(2.12519359588623) q[2];
sx q[2];
rz(4.3805526812845) q[2];
sx q[2];
rz(8.48003551959201) q[2];
rz(0.384812086820602) q[3];
sx q[3];
rz(2.02315691311891) q[3];
sx q[3];
rz(8.46772131919071) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.621895730495453) q[0];
sx q[0];
rz(5.21749392350251) q[0];
sx q[0];
rz(10.9863952159803) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-0.524295389652252) q[1];
sx q[1];
rz(1.36846628983552) q[1];
sx q[1];
rz(10.1069585442464) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-2.20484590530396) q[2];
sx q[2];
rz(4.31631949742372) q[2];
sx q[2];
rz(11.8699183225553) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(2.11123585700989) q[3];
sx q[3];
rz(4.64666417439515) q[3];
sx q[3];
rz(8.32820079325839) q[3];
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