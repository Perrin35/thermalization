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
rz(-0.090341180562973) q[0];
sx q[0];
rz(2.098896535235) q[0];
sx q[0];
rz(12.2876956224362) q[0];
rz(4.28840446472168) q[1];
sx q[1];
rz(0.524622591333934) q[1];
sx q[1];
rz(8.15370318888828) q[1];
cx q[1],q[0];
rz(-4.02567386627197) q[0];
sx q[0];
rz(1.83130982716615) q[0];
sx q[0];
rz(9.44367284550473) q[0];
rz(3.99568462371826) q[2];
sx q[2];
rz(6.08459726174409) q[2];
sx q[2];
rz(9.70303092002078) q[2];
cx q[2],q[1];
rz(-0.482857406139374) q[1];
sx q[1];
rz(7.00037685235078) q[1];
sx q[1];
rz(15.5871591329496) q[1];
rz(0.579234898090363) q[3];
sx q[3];
rz(2.0982200225168) q[3];
sx q[3];
rz(9.61949746905967) q[3];
cx q[3],q[2];
rz(-3.03113865852356) q[2];
sx q[2];
rz(5.06567827065522) q[2];
sx q[2];
rz(12.3041646242063) q[2];
rz(1.05552220344543) q[3];
sx q[3];
rz(4.388975532847) q[3];
sx q[3];
rz(9.51211566328212) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.5259040594101) q[0];
sx q[0];
rz(3.64510509570176) q[0];
sx q[0];
rz(9.2431887447755) q[0];
rz(4.88237476348877) q[1];
sx q[1];
rz(4.33435777028138) q[1];
sx q[1];
rz(8.58385852574512) q[1];
cx q[1],q[0];
rz(-1.4940482378006) q[0];
sx q[0];
rz(7.53398671944673) q[0];
sx q[0];
rz(12.5849432706754) q[0];
rz(-3.78149628639221) q[2];
sx q[2];
rz(8.24555030663545) q[2];
sx q[2];
rz(10.6118572711866) q[2];
cx q[2],q[1];
rz(-1.72791981697083) q[1];
sx q[1];
rz(0.172197254496165) q[1];
sx q[1];
rz(9.21212849616214) q[1];
rz(0.612003028392792) q[3];
sx q[3];
rz(5.01796940167482) q[3];
sx q[3];
rz(15.2540759801786) q[3];
cx q[3],q[2];
rz(5.90095663070679) q[2];
sx q[2];
rz(0.214919241266795) q[2];
sx q[2];
rz(14.6430735349576) q[2];
rz(6.22175645828247) q[3];
sx q[3];
rz(1.80627849896485) q[3];
sx q[3];
rz(11.0646841287534) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.43997287750244) q[0];
sx q[0];
rz(1.96236625512178) q[0];
sx q[0];
rz(9.22572060524627) q[0];
rz(0.82636034488678) q[1];
sx q[1];
rz(7.20143237908418) q[1];
sx q[1];
rz(7.0602130651395) q[1];
cx q[1],q[0];
rz(-0.137644052505493) q[0];
sx q[0];
rz(3.31193347473676) q[0];
sx q[0];
rz(6.68385932444736) q[0];
rz(-1.14116823673248) q[2];
sx q[2];
rz(5.21504607995088) q[2];
sx q[2];
rz(6.99188826083347) q[2];
cx q[2],q[1];
rz(-1.22328448295593) q[1];
sx q[1];
rz(4.45377698739106) q[1];
sx q[1];
rz(11.9997730016629) q[1];
rz(-0.593624591827393) q[3];
sx q[3];
rz(1.80261734326417) q[3];
sx q[3];
rz(7.32340476512119) q[3];
cx q[3],q[2];
rz(-1.10211539268494) q[2];
sx q[2];
rz(1.14810434182221) q[2];
sx q[2];
rz(12.0533284902494) q[2];
rz(-2.3593864440918) q[3];
sx q[3];
rz(4.34738090832765) q[3];
sx q[3];
rz(7.66006968020602) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.9272598028183) q[0];
sx q[0];
rz(4.73972037632997) q[0];
sx q[0];
rz(7.77367446421787) q[0];
rz(2.60930418968201) q[1];
sx q[1];
rz(2.51063010294969) q[1];
sx q[1];
rz(10.9623342513959) q[1];
cx q[1],q[0];
rz(-0.101899422705173) q[0];
sx q[0];
rz(2.94657747645909) q[0];
sx q[0];
rz(7.27443072795078) q[0];
rz(0.152875259518623) q[2];
sx q[2];
rz(5.73174563248689) q[2];
sx q[2];
rz(8.07514879702731) q[2];
cx q[2],q[1];
rz(3.72663736343384) q[1];
sx q[1];
rz(5.34661331971223) q[1];
sx q[1];
rz(10.9076246976773) q[1];
rz(0.234611481428146) q[3];
sx q[3];
rz(5.22777596314485) q[3];
sx q[3];
rz(8.38083682059451) q[3];
cx q[3],q[2];
rz(-1.60604667663574) q[2];
sx q[2];
rz(5.73637762864167) q[2];
sx q[2];
rz(9.43766302093073) q[2];
rz(-1.20312011241913) q[3];
sx q[3];
rz(4.60834530194337) q[3];
sx q[3];
rz(11.5446460008542) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.13836979866028) q[0];
sx q[0];
rz(2.51651814778382) q[0];
sx q[0];
rz(11.1665776729505) q[0];
rz(2.78955817222595) q[1];
sx q[1];
rz(2.08759179909761) q[1];
sx q[1];
rz(8.58693686722919) q[1];
cx q[1],q[0];
rz(-0.452338814735413) q[0];
sx q[0];
rz(4.86429384549195) q[0];
sx q[0];
rz(11.8041939496915) q[0];
rz(2.18508195877075) q[2];
sx q[2];
rz(2.09622422059114) q[2];
sx q[2];
rz(12.5050039052884) q[2];
cx q[2],q[1];
rz(-1.69882357120514) q[1];
sx q[1];
rz(4.83923724492128) q[1];
sx q[1];
rz(11.5640570878904) q[1];
rz(-0.153763055801392) q[3];
sx q[3];
rz(4.54991784890229) q[3];
sx q[3];
rz(10.211639738075) q[3];
cx q[3],q[2];
rz(1.02485227584839) q[2];
sx q[2];
rz(4.10631802876527) q[2];
sx q[2];
rz(10.6084508657376) q[2];
rz(0.270703256130219) q[3];
sx q[3];
rz(5.34271469910676) q[3];
sx q[3];
rz(11.0337378740232) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0266413316130638) q[0];
sx q[0];
rz(2.5178434570604) q[0];
sx q[0];
rz(11.9872295617978) q[0];
rz(-1.22224140167236) q[1];
sx q[1];
rz(1.83415833313996) q[1];
sx q[1];
rz(11.4567448854367) q[1];
cx q[1],q[0];
rz(0.189411759376526) q[0];
sx q[0];
rz(3.36643316050107) q[0];
sx q[0];
rz(10.2525615453641) q[0];
rz(1.26583981513977) q[2];
sx q[2];
rz(4.76622453530366) q[2];
sx q[2];
rz(10.1801560878675) q[2];
cx q[2],q[1];
rz(-3.06409692764282) q[1];
sx q[1];
rz(7.5228683074289) q[1];
sx q[1];
rz(12.8287536859433) q[1];
rz(1.05433130264282) q[3];
sx q[3];
rz(5.19599357445771) q[3];
sx q[3];
rz(11.6524019002835) q[3];
cx q[3],q[2];
rz(1.19607245922089) q[2];
sx q[2];
rz(4.98548713524873) q[2];
sx q[2];
rz(9.78223160504504) q[2];
rz(1.18999052047729) q[3];
sx q[3];
rz(4.34171214898164) q[3];
sx q[3];
rz(10.7644996404569) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.599151968955994) q[0];
sx q[0];
rz(4.15867200692231) q[0];
sx q[0];
rz(10.2168594360273) q[0];
rz(0.704108595848083) q[1];
sx q[1];
rz(6.73902145226533) q[1];
sx q[1];
rz(10.9832775354306) q[1];
cx q[1],q[0];
rz(1.96071994304657) q[0];
sx q[0];
rz(4.70944193203981) q[0];
sx q[0];
rz(11.4597642183225) q[0];
rz(0.75419819355011) q[2];
sx q[2];
rz(3.92023316224153) q[2];
sx q[2];
rz(13.739555335037) q[2];
cx q[2],q[1];
rz(1.66867911815643) q[1];
sx q[1];
rz(1.9967419226938) q[1];
sx q[1];
rz(9.38692688792154) q[1];
rz(2.32220268249512) q[3];
sx q[3];
rz(4.54611674149568) q[3];
sx q[3];
rz(10.0653224348943) q[3];
cx q[3],q[2];
rz(3.03029704093933) q[2];
sx q[2];
rz(0.788124950724193) q[2];
sx q[2];
rz(9.42896383245617) q[2];
rz(-0.266751021146774) q[3];
sx q[3];
rz(4.66560009320314) q[3];
sx q[3];
rz(8.70254693030521) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.580771446228027) q[0];
sx q[0];
rz(5.55367174943025) q[0];
sx q[0];
rz(9.58082289098903) q[0];
rz(-1.58496284484863) q[1];
sx q[1];
rz(2.27923509676988) q[1];
sx q[1];
rz(11.9716045617978) q[1];
cx q[1],q[0];
rz(-2.13463354110718) q[0];
sx q[0];
rz(1.31583205063874) q[0];
sx q[0];
rz(10.4862087726514) q[0];
rz(5.40472173690796) q[2];
sx q[2];
rz(4.93823734124238) q[2];
sx q[2];
rz(11.8747479677121) q[2];
cx q[2],q[1];
rz(-0.225603550672531) q[1];
sx q[1];
rz(4.02222404082353) q[1];
sx q[1];
rz(8.62119266986057) q[1];
rz(-0.849496245384216) q[3];
sx q[3];
rz(0.870660456018992) q[3];
sx q[3];
rz(6.96699998377963) q[3];
cx q[3],q[2];
rz(0.567383706569672) q[2];
sx q[2];
rz(4.70417872269685) q[2];
sx q[2];
rz(10.4527700900952) q[2];
rz(-1.79936408996582) q[3];
sx q[3];
rz(2.48643109400804) q[3];
sx q[3];
rz(10.7596013307492) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.87803554534912) q[0];
sx q[0];
rz(4.5981309731775) q[0];
sx q[0];
rz(8.626517212383) q[0];
rz(2.34106659889221) q[1];
sx q[1];
rz(3.63056749303872) q[1];
sx q[1];
rz(8.87854561804935) q[1];
cx q[1],q[0];
rz(3.67085886001587) q[0];
sx q[0];
rz(3.7533427794748) q[0];
sx q[0];
rz(12.3762111425321) q[0];
rz(0.428225487470627) q[2];
sx q[2];
rz(5.25830498536164) q[2];
sx q[2];
rz(8.52984706162616) q[2];
cx q[2],q[1];
rz(-0.479104667901993) q[1];
sx q[1];
rz(2.42913851340348) q[1];
sx q[1];
rz(11.7893621683042) q[1];
rz(1.26359748840332) q[3];
sx q[3];
rz(3.37874752481515) q[3];
sx q[3];
rz(12.4568550348203) q[3];
cx q[3],q[2];
rz(-0.747867822647095) q[2];
sx q[2];
rz(4.23168245156343) q[2];
sx q[2];
rz(9.36174292712613) q[2];
rz(0.697832703590393) q[3];
sx q[3];
rz(6.02501860459382) q[3];
sx q[3];
rz(8.51168713568851) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.495527118444443) q[0];
sx q[0];
rz(5.40554061730439) q[0];
sx q[0];
rz(9.00873596071407) q[0];
rz(1.3620194196701) q[1];
sx q[1];
rz(5.15905443032319) q[1];
sx q[1];
rz(10.7175664663236) q[1];
cx q[1],q[0];
rz(1.10716152191162) q[0];
sx q[0];
rz(3.49784809549386) q[0];
sx q[0];
rz(10.2866885423581) q[0];
rz(1.1146194934845) q[2];
sx q[2];
rz(4.96451547940309) q[2];
sx q[2];
rz(13.7415270566861) q[2];
cx q[2],q[1];
rz(2.46830081939697) q[1];
sx q[1];
rz(2.39299622376496) q[1];
sx q[1];
rz(10.1382016897123) q[1];
rz(-0.536968886852264) q[3];
sx q[3];
rz(6.76467314560945) q[3];
sx q[3];
rz(9.43531513436838) q[3];
cx q[3],q[2];
rz(-3.82459139823914) q[2];
sx q[2];
rz(4.47799709637696) q[2];
sx q[2];
rz(10.3322953343312) q[2];
rz(-1.26097428798676) q[3];
sx q[3];
rz(3.71615538199479) q[3];
sx q[3];
rz(11.1169615745465) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.818953454494476) q[0];
sx q[0];
rz(-0.837539521855764) q[0];
sx q[0];
rz(9.26753327845737) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-3.75367283821106) q[1];
sx q[1];
rz(5.83540788491304) q[1];
sx q[1];
rz(9.98050490616962) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.914084672927856) q[2];
sx q[2];
rz(5.16633907158906) q[2];
sx q[2];
rz(7.82577631472751) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.000672237714752555) q[3];
sx q[3];
rz(3.00816248555715) q[3];
sx q[3];
rz(7.12933990954562) q[3];
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
