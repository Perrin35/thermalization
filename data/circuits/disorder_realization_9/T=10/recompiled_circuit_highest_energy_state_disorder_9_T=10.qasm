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
rz(0.23586805164814) q[0];
sx q[0];
rz(1.996082456904) q[0];
sx q[0];
rz(9.37620147540375) q[0];
rz(1.51611685752869) q[1];
sx q[1];
rz(4.23443439801271) q[1];
sx q[1];
rz(8.70901272296115) q[1];
cx q[1],q[0];
rz(-0.193120345473289) q[0];
sx q[0];
rz(3.31491022010381) q[0];
sx q[0];
rz(9.59620348214313) q[0];
rz(-0.279996544122696) q[2];
sx q[2];
rz(4.91152337391908) q[2];
sx q[2];
rz(9.37382016926214) q[2];
cx q[2],q[1];
rz(-0.298083335161209) q[1];
sx q[1];
rz(4.02710268099839) q[1];
sx q[1];
rz(9.95047167538806) q[1];
rz(4.06800651550293) q[3];
sx q[3];
rz(3.75999936659867) q[3];
sx q[3];
rz(8.23909387587711) q[3];
cx q[3],q[2];
rz(0.488926589488983) q[2];
sx q[2];
rz(2.58161398966844) q[2];
sx q[2];
rz(10.4803114890973) q[2];
rz(-3.54759550094604) q[3];
sx q[3];
rz(4.45522001584107) q[3];
sx q[3];
rz(10.2463794708173) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.157762885093689) q[0];
sx q[0];
rz(2.10217479069764) q[0];
sx q[0];
rz(9.98095343112155) q[0];
rz(-1.32933855056763) q[1];
sx q[1];
rz(3.76032564242417) q[1];
sx q[1];
rz(12.4830975294034) q[1];
cx q[1],q[0];
rz(0.525169253349304) q[0];
sx q[0];
rz(0.942481907206126) q[0];
sx q[0];
rz(9.82798371314212) q[0];
rz(0.742342889308929) q[2];
sx q[2];
rz(3.90861681302125) q[2];
sx q[2];
rz(7.67567155360385) q[2];
cx q[2],q[1];
rz(-1.14948070049286) q[1];
sx q[1];
rz(3.64210298855836) q[1];
sx q[1];
rz(13.6652941465299) q[1];
rz(1.07636916637421) q[3];
sx q[3];
rz(5.22321739991242) q[3];
sx q[3];
rz(8.31470987795993) q[3];
cx q[3],q[2];
rz(-0.197815299034119) q[2];
sx q[2];
rz(1.75217202504212) q[2];
sx q[2];
rz(9.91754145025417) q[2];
rz(-1.02657759189606) q[3];
sx q[3];
rz(2.83870983322198) q[3];
sx q[3];
rz(11.153806424133) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.763484120368958) q[0];
sx q[0];
rz(3.70430258114869) q[0];
sx q[0];
rz(8.98599836825534) q[0];
rz(1.45250606536865) q[1];
sx q[1];
rz(5.96135249932344) q[1];
sx q[1];
rz(12.1734173059384) q[1];
cx q[1],q[0];
rz(1.90714359283447) q[0];
sx q[0];
rz(3.63357815344865) q[0];
sx q[0];
rz(9.48391114770576) q[0];
rz(3.19767189025879) q[2];
sx q[2];
rz(2.61255565484101) q[2];
sx q[2];
rz(9.03533390759631) q[2];
cx q[2],q[1];
rz(0.716869831085205) q[1];
sx q[1];
rz(3.56574392517144) q[1];
sx q[1];
rz(9.59798427521392) q[1];
rz(0.921513795852661) q[3];
sx q[3];
rz(4.77577212651307) q[3];
sx q[3];
rz(11.0646083116452) q[3];
cx q[3],q[2];
rz(0.689859330654144) q[2];
sx q[2];
rz(5.43875017960603) q[2];
sx q[2];
rz(9.15125099419757) q[2];
rz(0.612559735774994) q[3];
sx q[3];
rz(1.72756830056245) q[3];
sx q[3];
rz(10.2812684535901) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.813711524009705) q[0];
sx q[0];
rz(3.73469081719453) q[0];
sx q[0];
rz(10.2255284547727) q[0];
rz(0.712093591690063) q[1];
sx q[1];
rz(4.26158228714997) q[1];
sx q[1];
rz(9.90794274806186) q[1];
cx q[1],q[0];
rz(1.90948510169983) q[0];
sx q[0];
rz(4.71299520333345) q[0];
sx q[0];
rz(8.59667954444095) q[0];
rz(2.64473414421082) q[2];
sx q[2];
rz(4.2261005957895) q[2];
sx q[2];
rz(8.69502005576297) q[2];
cx q[2],q[1];
rz(2.93127155303955) q[1];
sx q[1];
rz(0.646697910624095) q[1];
sx q[1];
rz(7.90184209345981) q[1];
rz(-0.336709529161453) q[3];
sx q[3];
rz(4.05129882891709) q[3];
sx q[3];
rz(9.11236894725963) q[3];
cx q[3],q[2];
rz(0.983727693557739) q[2];
sx q[2];
rz(3.85220769246156) q[2];
sx q[2];
rz(10.7092896461408) q[2];
rz(0.487867832183838) q[3];
sx q[3];
rz(4.84593871434266) q[3];
sx q[3];
rz(10.8792639732282) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.812164127826691) q[0];
sx q[0];
rz(2.28211304743821) q[0];
sx q[0];
rz(12.101722931854) q[0];
rz(5.33477115631104) q[1];
sx q[1];
rz(2.30348381598527) q[1];
sx q[1];
rz(7.9714040517728) q[1];
cx q[1],q[0];
rz(-0.689357161521912) q[0];
sx q[0];
rz(3.18472980906303) q[0];
sx q[0];
rz(8.79568318127795) q[0];
rz(0.830346882343292) q[2];
sx q[2];
rz(4.5210749228769) q[2];
sx q[2];
rz(11.6054389238279) q[2];
cx q[2],q[1];
rz(1.59111154079437) q[1];
sx q[1];
rz(4.77948621113832) q[1];
sx q[1];
rz(10.3842068672101) q[1];
rz(2.29591345787048) q[3];
sx q[3];
rz(4.31021025975282) q[3];
sx q[3];
rz(8.45652202366992) q[3];
cx q[3],q[2];
rz(-0.984430849552155) q[2];
sx q[2];
rz(4.18899455864961) q[2];
sx q[2];
rz(9.0303380548875) q[2];
rz(-1.98583424091339) q[3];
sx q[3];
rz(4.07164934475953) q[3];
sx q[3];
rz(11.2578927040021) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.30094575881958) q[0];
sx q[0];
rz(3.55738946993882) q[0];
sx q[0];
rz(10.4867104053418) q[0];
rz(-0.951470792293549) q[1];
sx q[1];
rz(4.01522609789903) q[1];
sx q[1];
rz(10.945509171478) q[1];
cx q[1],q[0];
rz(0.0117640476673841) q[0];
sx q[0];
rz(4.7655589898401) q[0];
sx q[0];
rz(12.2891421079557) q[0];
rz(0.925918400287628) q[2];
sx q[2];
rz(5.19038501580293) q[2];
sx q[2];
rz(8.42419645785495) q[2];
cx q[2],q[1];
rz(0.86506062746048) q[1];
sx q[1];
rz(3.92621961434419) q[1];
sx q[1];
rz(8.08577618598148) q[1];
rz(2.06534814834595) q[3];
sx q[3];
rz(2.7003717144304) q[3];
sx q[3];
rz(8.38210568427249) q[3];
cx q[3],q[2];
rz(1.11573970317841) q[2];
sx q[2];
rz(1.91511765320832) q[2];
sx q[2];
rz(9.17907304166957) q[2];
rz(-1.45415163040161) q[3];
sx q[3];
rz(4.77155569394166) q[3];
sx q[3];
rz(10.2551264524381) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.655697345733643) q[0];
sx q[0];
rz(2.40351543028886) q[0];
sx q[0];
rz(10.1194990634839) q[0];
rz(3.04095363616943) q[1];
sx q[1];
rz(5.24742451508576) q[1];
sx q[1];
rz(7.14839384555026) q[1];
cx q[1],q[0];
rz(0.145102053880692) q[0];
sx q[0];
rz(4.29443398316438) q[0];
sx q[0];
rz(9.47017187102839) q[0];
rz(0.660431921482086) q[2];
sx q[2];
rz(4.72730604012544) q[2];
sx q[2];
rz(7.21494553088352) q[2];
cx q[2],q[1];
rz(0.455649018287659) q[1];
sx q[1];
rz(5.65075436432893) q[1];
sx q[1];
rz(9.16976380943462) q[1];
rz(0.551366150379181) q[3];
sx q[3];
rz(2.16882810195024) q[3];
sx q[3];
rz(10.1674294233243) q[3];
cx q[3],q[2];
rz(0.184090986847878) q[2];
sx q[2];
rz(3.89506599505479) q[2];
sx q[2];
rz(8.51457319258853) q[2];
rz(1.45220673084259) q[3];
sx q[3];
rz(1.9096749146753) q[3];
sx q[3];
rz(8.96049544810458) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.243773952126503) q[0];
sx q[0];
rz(2.09385660489137) q[0];
sx q[0];
rz(8.80149564742252) q[0];
rz(-0.155776232481003) q[1];
sx q[1];
rz(3.91889026959474) q[1];
sx q[1];
rz(9.16146937607929) q[1];
cx q[1],q[0];
rz(2.14151883125305) q[0];
sx q[0];
rz(2.92630683084066) q[0];
sx q[0];
rz(9.3747100405316) q[0];
rz(0.901562988758087) q[2];
sx q[2];
rz(3.04416360159452) q[2];
sx q[2];
rz(11.6395239591519) q[2];
cx q[2],q[1];
rz(4.44310331344604) q[1];
sx q[1];
rz(3.7231908758455) q[1];
sx q[1];
rz(7.71368870734378) q[1];
rz(0.543680012226105) q[3];
sx q[3];
rz(4.00461349089677) q[3];
sx q[3];
rz(10.3100279331128) q[3];
cx q[3],q[2];
rz(0.891736328601837) q[2];
sx q[2];
rz(2.89204341371591) q[2];
sx q[2];
rz(10.3128311991613) q[2];
rz(2.21372103691101) q[3];
sx q[3];
rz(4.26952317555482) q[3];
sx q[3];
rz(8.45169148444339) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.718582808971405) q[0];
sx q[0];
rz(3.65123745997483) q[0];
sx q[0];
rz(8.92404124735996) q[0];
rz(2.02107334136963) q[1];
sx q[1];
rz(5.50214591820771) q[1];
sx q[1];
rz(10.226263499252) q[1];
cx q[1],q[0];
rz(0.0285417381674051) q[0];
sx q[0];
rz(3.07059122820432) q[0];
sx q[0];
rz(7.82619116305515) q[0];
rz(-0.525252103805542) q[2];
sx q[2];
rz(4.21039024193818) q[2];
sx q[2];
rz(10.5096401929776) q[2];
cx q[2],q[1];
rz(1.55133306980133) q[1];
sx q[1];
rz(4.09671095212037) q[1];
sx q[1];
rz(8.89214465617343) q[1];
rz(0.257973998785019) q[3];
sx q[3];
rz(4.52387908299501) q[3];
sx q[3];
rz(9.21352080105945) q[3];
cx q[3],q[2];
rz(0.720061719417572) q[2];
sx q[2];
rz(4.46300330956513) q[2];
sx q[2];
rz(12.3361751794736) q[2];
rz(-0.00975333340466022) q[3];
sx q[3];
rz(4.7662582715326) q[3];
sx q[3];
rz(9.58074176906749) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.44740509986877) q[0];
sx q[0];
rz(4.91560128529603) q[0];
sx q[0];
rz(10.194378232948) q[0];
rz(-2.17904782295227) q[1];
sx q[1];
rz(4.45943346818025) q[1];
sx q[1];
rz(11.5792751073758) q[1];
cx q[1],q[0];
rz(2.01131629943848) q[0];
sx q[0];
rz(3.67302218277986) q[0];
sx q[0];
rz(9.56269606053039) q[0];
rz(-0.55743420124054) q[2];
sx q[2];
rz(3.42562806804711) q[2];
sx q[2];
rz(10.782383298866) q[2];
cx q[2],q[1];
rz(-1.92224395275116) q[1];
sx q[1];
rz(4.74551394780213) q[1];
sx q[1];
rz(9.38354404865905) q[1];
rz(-0.494786649942398) q[3];
sx q[3];
rz(5.42519393761689) q[3];
sx q[3];
rz(9.53515376000806) q[3];
cx q[3],q[2];
rz(1.71396815776825) q[2];
sx q[2];
rz(4.07887402375276) q[2];
sx q[2];
rz(9.49289981125995) q[2];
rz(-0.449759751558304) q[3];
sx q[3];
rz(3.55683049758012) q[3];
sx q[3];
rz(10.8629818916242) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.208520516753197) q[0];
sx q[0];
rz(3.31934896309907) q[0];
sx q[0];
rz(9.64129520057842) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(1.76726794242859) q[1];
sx q[1];
rz(4.92294541199739) q[1];
sx q[1];
rz(10.6192935466687) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(2.12404870986938) q[2];
sx q[2];
rz(4.44718781312043) q[2];
sx q[2];
rz(7.33079192637607) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.979496657848358) q[3];
sx q[3];
rz(2.50801620085771) q[3];
sx q[3];
rz(10.1393651723783) q[3];
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
