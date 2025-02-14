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
rz(0.718304634094238) q[0];
sx q[0];
rz(2.36513200600679) q[0];
sx q[0];
rz(9.79070103763744) q[0];
rz(-2.79468202590942) q[1];
sx q[1];
rz(3.42630979617173) q[1];
sx q[1];
rz(10.8135362625043) q[1];
cx q[1],q[0];
rz(0.589088678359985) q[0];
sx q[0];
rz(3.35461071332032) q[0];
sx q[0];
rz(11.3159826755445) q[0];
rz(2.26348447799683) q[2];
sx q[2];
rz(4.48035493691499) q[2];
sx q[2];
rz(7.50499186515018) q[2];
cx q[2],q[1];
rz(1.54688155651093) q[1];
sx q[1];
rz(4.32339456875856) q[1];
sx q[1];
rz(8.94238457678958) q[1];
rz(3.99929523468018) q[3];
sx q[3];
rz(3.8700746019655) q[3];
sx q[3];
rz(8.4301892876546) q[3];
cx q[3],q[2];
rz(-2.18467879295349) q[2];
sx q[2];
rz(4.06173667510087) q[2];
sx q[2];
rz(13.5015682935636) q[2];
rz(-0.301855593919754) q[3];
sx q[3];
rz(4.7409018595987) q[3];
sx q[3];
rz(9.03111711739703) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.16441559791565) q[0];
sx q[0];
rz(1.26576069195802) q[0];
sx q[0];
rz(8.92648301123782) q[0];
rz(-1.71388185024261) q[1];
sx q[1];
rz(1.93528130848939) q[1];
sx q[1];
rz(10.5115364551465) q[1];
cx q[1],q[0];
rz(0.552884936332703) q[0];
sx q[0];
rz(3.87947145302827) q[0];
sx q[0];
rz(14.4928121328275) q[0];
rz(-0.266579657793045) q[2];
sx q[2];
rz(4.89949563344056) q[2];
sx q[2];
rz(9.53026705085441) q[2];
cx q[2],q[1];
rz(5.81049346923828) q[1];
sx q[1];
rz(7.69367185433442) q[1];
sx q[1];
rz(7.41326782702609) q[1];
rz(-1.60727536678314) q[3];
sx q[3];
rz(3.46000561316545) q[3];
sx q[3];
rz(8.48088542222186) q[3];
cx q[3],q[2];
rz(3.34024357795715) q[2];
sx q[2];
rz(5.9725765307718) q[2];
sx q[2];
rz(14.3757419347684) q[2];
rz(-1.99409353733063) q[3];
sx q[3];
rz(4.70445266564424) q[3];
sx q[3];
rz(8.09214804171726) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.93208169937134) q[0];
sx q[0];
rz(1.40283432801301) q[0];
sx q[0];
rz(6.12594864367648) q[0];
rz(1.22781848907471) q[1];
sx q[1];
rz(2.93240517576272) q[1];
sx q[1];
rz(9.00270617603465) q[1];
cx q[1],q[0];
rz(-1.87647187709808) q[0];
sx q[0];
rz(2.94199638267095) q[0];
sx q[0];
rz(8.95467964409992) q[0];
rz(-1.29319620132446) q[2];
sx q[2];
rz(4.54500082333619) q[2];
sx q[2];
rz(12.5244278669278) q[2];
cx q[2],q[1];
rz(-1.28733086585999) q[1];
sx q[1];
rz(3.53347701032693) q[1];
sx q[1];
rz(12.3336839437406) q[1];
rz(-3.45895624160767) q[3];
sx q[3];
rz(3.93606397707994) q[3];
sx q[3];
rz(10.4251850604932) q[3];
cx q[3],q[2];
rz(-0.328553646802902) q[2];
sx q[2];
rz(5.32287207444245) q[2];
sx q[2];
rz(10.1588462352674) q[2];
rz(1.07753932476044) q[3];
sx q[3];
rz(2.45340928633744) q[3];
sx q[3];
rz(12.49116060733) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.58287358283997) q[0];
sx q[0];
rz(4.15210834343965) q[0];
sx q[0];
rz(6.29966661929294) q[0];
rz(3.21613550186157) q[1];
sx q[1];
rz(2.68389791448648) q[1];
sx q[1];
rz(9.16963887809917) q[1];
cx q[1],q[0];
rz(-0.821273863315582) q[0];
sx q[0];
rz(3.6458285172754) q[0];
sx q[0];
rz(8.1551244020383) q[0];
rz(1.90569126605988) q[2];
sx q[2];
rz(4.17195585568482) q[2];
sx q[2];
rz(10.6135139226834) q[2];
cx q[2],q[1];
rz(-0.296834766864777) q[1];
sx q[1];
rz(2.18983087142045) q[1];
sx q[1];
rz(12.0311882257383) q[1];
rz(1.68473720550537) q[3];
sx q[3];
rz(4.36503210862214) q[3];
sx q[3];
rz(7.52641234397098) q[3];
cx q[3],q[2];
rz(-2.52693462371826) q[2];
sx q[2];
rz(5.58607545693452) q[2];
sx q[2];
rz(10.1245796441953) q[2];
rz(-3.049729347229) q[3];
sx q[3];
rz(5.06586876709992) q[3];
sx q[3];
rz(9.87950382231876) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.622345387935638) q[0];
sx q[0];
rz(5.46280279954011) q[0];
sx q[0];
rz(10.5366111755292) q[0];
rz(0.190192922949791) q[1];
sx q[1];
rz(5.39804282982881) q[1];
sx q[1];
rz(10.5846518039624) q[1];
cx q[1],q[0];
rz(-0.969441890716553) q[0];
sx q[0];
rz(4.74686518509919) q[0];
sx q[0];
rz(12.9102496862332) q[0];
rz(0.196148857474327) q[2];
sx q[2];
rz(3.7892394383722) q[2];
sx q[2];
rz(12.1886212587278) q[2];
cx q[2],q[1];
rz(-0.240073561668396) q[1];
sx q[1];
rz(4.1699881871515) q[1];
sx q[1];
rz(12.276066517822) q[1];
rz(1.61574792861938) q[3];
sx q[3];
rz(4.19093564351136) q[3];
sx q[3];
rz(11.002766108505) q[3];
cx q[3],q[2];
rz(1.57088923454285) q[2];
sx q[2];
rz(5.46001520951326) q[2];
sx q[2];
rz(9.19432043134376) q[2];
rz(1.06200909614563) q[3];
sx q[3];
rz(4.41746500332887) q[3];
sx q[3];
rz(7.91637108325168) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.452881693840027) q[0];
sx q[0];
rz(2.08275607426698) q[0];
sx q[0];
rz(11.2105091571729) q[0];
rz(2.61182498931885) q[1];
sx q[1];
rz(3.92387655575807) q[1];
sx q[1];
rz(11.4145381212155) q[1];
cx q[1],q[0];
rz(0.0655819475650787) q[0];
sx q[0];
rz(4.35706499417359) q[0];
sx q[0];
rz(9.53003978579446) q[0];
rz(0.157635018229485) q[2];
sx q[2];
rz(5.73602047761018) q[2];
sx q[2];
rz(13.1947982072751) q[2];
cx q[2],q[1];
rz(-0.639101445674896) q[1];
sx q[1];
rz(4.51848998864228) q[1];
sx q[1];
rz(10.6246049165647) q[1];
rz(-0.117529585957527) q[3];
sx q[3];
rz(4.4995974620157) q[3];
sx q[3];
rz(12.8301045656125) q[3];
cx q[3],q[2];
rz(-1.66662240028381) q[2];
sx q[2];
rz(5.64184871514375) q[2];
sx q[2];
rz(8.95767638682529) q[2];
rz(1.01116716861725) q[3];
sx q[3];
rz(5.93953147728974) q[3];
sx q[3];
rz(10.7941971778791) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.29246270656586) q[0];
sx q[0];
rz(0.155015381174632) q[0];
sx q[0];
rz(12.3416840791623) q[0];
rz(-3.30541229248047) q[1];
sx q[1];
rz(5.94406357605989) q[1];
sx q[1];
rz(6.92954108714267) q[1];
cx q[1],q[0];
rz(-0.408989727497101) q[0];
sx q[0];
rz(1.32091310818727) q[0];
sx q[0];
rz(8.88215336798831) q[0];
rz(0.322856307029724) q[2];
sx q[2];
rz(4.50528553326661) q[2];
sx q[2];
rz(9.75237975119754) q[2];
cx q[2],q[1];
rz(0.0283706281334162) q[1];
sx q[1];
rz(5.73444763024385) q[1];
sx q[1];
rz(10.8509065866391) q[1];
rz(1.89432370662689) q[3];
sx q[3];
rz(3.54801112611825) q[3];
sx q[3];
rz(12.9930002450864) q[3];
cx q[3],q[2];
rz(-0.155676379799843) q[2];
sx q[2];
rz(4.83515837987001) q[2];
sx q[2];
rz(9.83473405837222) q[2];
rz(3.9719386100769) q[3];
sx q[3];
rz(4.09032342036302) q[3];
sx q[3];
rz(7.73108396529361) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.537725627422333) q[0];
sx q[0];
rz(5.58997646172578) q[0];
sx q[0];
rz(9.17847334443733) q[0];
rz(3.96818828582764) q[1];
sx q[1];
rz(1.39845028718049) q[1];
sx q[1];
rz(9.68873891829654) q[1];
cx q[1],q[0];
rz(1.95190989971161) q[0];
sx q[0];
rz(6.37893215020234) q[0];
sx q[0];
rz(9.34907632171317) q[0];
rz(2.3545458316803) q[2];
sx q[2];
rz(5.06199851830537) q[2];
sx q[2];
rz(9.14905340074703) q[2];
cx q[2],q[1];
rz(-1.82496559619904) q[1];
sx q[1];
rz(2.39676031668717) q[1];
sx q[1];
rz(9.00321677922412) q[1];
rz(-0.983718693256378) q[3];
sx q[3];
rz(4.34285322030122) q[3];
sx q[3];
rz(7.33042738436862) q[3];
cx q[3],q[2];
rz(3.24653387069702) q[2];
sx q[2];
rz(1.87373259862001) q[2];
sx q[2];
rz(8.70966056584522) q[2];
rz(0.0713083297014236) q[3];
sx q[3];
rz(4.72619858582551) q[3];
sx q[3];
rz(7.030070757858) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.00953733921051) q[0];
sx q[0];
rz(1.65224638779695) q[0];
sx q[0];
rz(12.1308309793393) q[0];
rz(-0.147773310542107) q[1];
sx q[1];
rz(4.85158196290071) q[1];
sx q[1];
rz(10.8933423519055) q[1];
cx q[1],q[0];
rz(-0.763984739780426) q[0];
sx q[0];
rz(5.14335623581941) q[0];
sx q[0];
rz(9.62280999719306) q[0];
rz(1.66506171226501) q[2];
sx q[2];
rz(5.58361163933808) q[2];
sx q[2];
rz(9.88110617398425) q[2];
cx q[2],q[1];
rz(0.109941534698009) q[1];
sx q[1];
rz(5.09958973725373) q[1];
sx q[1];
rz(11.1757811069409) q[1];
rz(-1.09575593471527) q[3];
sx q[3];
rz(6.73092547257478) q[3];
sx q[3];
rz(10.489295577995) q[3];
cx q[3],q[2];
rz(-0.744641661643982) q[2];
sx q[2];
rz(8.04725948174531) q[2];
sx q[2];
rz(7.17612407206699) q[2];
rz(1.86301553249359) q[3];
sx q[3];
rz(5.1262927373224) q[3];
sx q[3];
rz(10.8881871461789) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.294991493225098) q[0];
sx q[0];
rz(5.95201936562593) q[0];
sx q[0];
rz(8.81871817111179) q[0];
rz(0.010295707732439) q[1];
sx q[1];
rz(2.1538486798578) q[1];
sx q[1];
rz(9.50470267831489) q[1];
cx q[1],q[0];
rz(2.85113596916199) q[0];
sx q[0];
rz(6.9093612750345) q[0];
sx q[0];
rz(9.86400390266582) q[0];
rz(2.41315793991089) q[2];
sx q[2];
rz(4.92951479752595) q[2];
sx q[2];
rz(12.1170022249143) q[2];
cx q[2],q[1];
rz(0.273516803979874) q[1];
sx q[1];
rz(3.54925257165963) q[1];
sx q[1];
rz(11.0516282081525) q[1];
rz(0.726945996284485) q[3];
sx q[3];
rz(5.97876396973664) q[3];
sx q[3];
rz(10.1812819600026) q[3];
cx q[3],q[2];
rz(3.00048780441284) q[2];
sx q[2];
rz(4.06844839652116) q[2];
sx q[2];
rz(7.29837915896579) q[2];
rz(0.601239323616028) q[3];
sx q[3];
rz(2.43981179793412) q[3];
sx q[3];
rz(11.5198220968167) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.970305025577545) q[0];
sx q[0];
rz(5.39873638947541) q[0];
sx q[0];
rz(9.49373484253093) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(2.52189230918884) q[1];
sx q[1];
rz(4.43505254586274) q[1];
sx q[1];
rz(8.48002103566333) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(2.82439064979553) q[2];
sx q[2];
rz(3.7404404600435) q[2];
sx q[2];
rz(6.12862203120395) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-1.06425309181213) q[3];
sx q[3];
rz(0.781305464106151) q[3];
sx q[3];
rz(11.2781688928525) q[3];
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
