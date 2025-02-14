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
rz(3.00388598442078) q[0];
sx q[0];
rz(0.51578155358369) q[0];
sx q[0];
rz(9.46573822050496) q[0];
rz(0.750866830348969) q[1];
sx q[1];
rz(2.65843957861001) q[1];
sx q[1];
rz(9.2106341779153) q[1];
cx q[1],q[0];
rz(1.0182603597641) q[0];
sx q[0];
rz(1.40544477303559) q[0];
sx q[0];
rz(9.33425257950231) q[0];
rz(-0.198039144277573) q[2];
sx q[2];
rz(1.65743401845033) q[2];
sx q[2];
rz(10.9660665750425) q[2];
cx q[2],q[1];
rz(-4.25351333618164) q[1];
sx q[1];
rz(2.40555337269837) q[1];
sx q[1];
rz(15.9379501104276) q[1];
rz(2.05841875076294) q[3];
sx q[3];
rz(2.02425256569917) q[3];
sx q[3];
rz(7.08886239527866) q[3];
cx q[3],q[2];
rz(0.21138770878315) q[2];
sx q[2];
rz(-1.38964208762114) q[2];
sx q[2];
rz(7.28156588076755) q[2];
rz(-1.80066478252411) q[3];
sx q[3];
rz(5.212626131373) q[3];
sx q[3];
rz(10.0827374815862) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.200876548886299) q[0];
sx q[0];
rz(4.05439606507356) q[0];
sx q[0];
rz(10.3024692892949) q[0];
rz(3.3974130153656) q[1];
sx q[1];
rz(4.10113039811189) q[1];
sx q[1];
rz(8.97604543565913) q[1];
cx q[1],q[0];
rz(3.43397045135498) q[0];
sx q[0];
rz(1.51499822934205) q[0];
sx q[0];
rz(10.6660154819409) q[0];
rz(-3.35541319847107) q[2];
sx q[2];
rz(4.81752112706239) q[2];
sx q[2];
rz(11.6563544034879) q[2];
cx q[2],q[1];
rz(-0.70905590057373) q[1];
sx q[1];
rz(4.01767978270585) q[1];
sx q[1];
rz(6.72880408763095) q[1];
rz(1.12345480918884) q[3];
sx q[3];
rz(5.29956379731233) q[3];
sx q[3];
rz(9.25472802519008) q[3];
cx q[3],q[2];
rz(-1.14827072620392) q[2];
sx q[2];
rz(1.69308701355989) q[2];
sx q[2];
rz(11.0235744476239) q[2];
rz(0.791050612926483) q[3];
sx q[3];
rz(1.4619210084253) q[3];
sx q[3];
rz(10.2171231865804) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.17447710037231) q[0];
sx q[0];
rz(3.73456904490525) q[0];
sx q[0];
rz(10.7052823066632) q[0];
rz(3.35060286521912) q[1];
sx q[1];
rz(5.1447321494394) q[1];
sx q[1];
rz(7.99762580393955) q[1];
cx q[1],q[0];
rz(-1.20365941524506) q[0];
sx q[0];
rz(5.26409211953218) q[0];
sx q[0];
rz(15.1791510343473) q[0];
rz(0.572978496551514) q[2];
sx q[2];
rz(6.48183050950105) q[2];
sx q[2];
rz(7.63261196612521) q[2];
cx q[2],q[1];
rz(1.18759655952454) q[1];
sx q[1];
rz(5.30203381379182) q[1];
sx q[1];
rz(9.7512436568658) q[1];
rz(-0.986737012863159) q[3];
sx q[3];
rz(4.82148793538148) q[3];
sx q[3];
rz(11.4857284784238) q[3];
cx q[3],q[2];
rz(-1.78724944591522) q[2];
sx q[2];
rz(2.30865827401216) q[2];
sx q[2];
rz(8.42616698741123) q[2];
rz(2.78658151626587) q[3];
sx q[3];
rz(4.52549269993836) q[3];
sx q[3];
rz(11.342083311073) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.760964751243591) q[0];
sx q[0];
rz(0.172920854883738) q[0];
sx q[0];
rz(10.3372689843099) q[0];
rz(-1.64187049865723) q[1];
sx q[1];
rz(4.38760975201661) q[1];
sx q[1];
rz(7.9610680103223) q[1];
cx q[1],q[0];
rz(2.64482235908508) q[0];
sx q[0];
rz(3.45120683510835) q[0];
sx q[0];
rz(13.792106127731) q[0];
rz(-3.05663752555847) q[2];
sx q[2];
rz(2.17329970200593) q[2];
sx q[2];
rz(9.2143497377555) q[2];
cx q[2],q[1];
rz(1.15965020656586) q[1];
sx q[1];
rz(4.4809434731775) q[1];
sx q[1];
rz(4.31369969844028) q[1];
rz(0.312954515218735) q[3];
sx q[3];
rz(3.88872906764085) q[3];
sx q[3];
rz(11.4299835920255) q[3];
cx q[3],q[2];
rz(4.203209400177) q[2];
sx q[2];
rz(1.43922630150849) q[2];
sx q[2];
rz(8.03269884585544) q[2];
rz(1.12130868434906) q[3];
sx q[3];
rz(1.08713546593721) q[3];
sx q[3];
rz(9.4794969022195) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0120809720829129) q[0];
sx q[0];
rz(4.71660724480683) q[0];
sx q[0];
rz(8.94681543707057) q[0];
rz(-1.67320084571838) q[1];
sx q[1];
rz(4.71123257477815) q[1];
sx q[1];
rz(10.4881089687268) q[1];
cx q[1],q[0];
rz(-0.376035183668137) q[0];
sx q[0];
rz(2.40046939452226) q[0];
sx q[0];
rz(11.8222560644071) q[0];
rz(9.188720703125) q[2];
sx q[2];
rz(6.85277620156343) q[2];
sx q[2];
rz(12.0131087064664) q[2];
cx q[2],q[1];
rz(-4.98363780975342) q[1];
sx q[1];
rz(5.40474501450593) q[1];
sx q[1];
rz(5.73015973567172) q[1];
rz(0.786126017570496) q[3];
sx q[3];
rz(5.86834731896455) q[3];
sx q[3];
rz(11.3945107221524) q[3];
cx q[3],q[2];
rz(0.310745060443878) q[2];
sx q[2];
rz(5.1060461123758) q[2];
sx q[2];
rz(11.3177184820096) q[2];
rz(2.25110745429993) q[3];
sx q[3];
rz(4.36314073403413) q[3];
sx q[3];
rz(9.85209388136073) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.86521995067596) q[0];
sx q[0];
rz(5.89738384087617) q[0];
sx q[0];
rz(10.6307877063672) q[0];
rz(1.68458664417267) q[1];
sx q[1];
rz(5.28866163094575) q[1];
sx q[1];
rz(10.438367342941) q[1];
cx q[1],q[0];
rz(0.591313660144806) q[0];
sx q[0];
rz(6.3489889224344) q[0];
sx q[0];
rz(3.99814746379062) q[0];
rz(3.57421064376831) q[2];
sx q[2];
rz(4.16574576695497) q[2];
sx q[2];
rz(8.34631762503787) q[2];
cx q[2],q[1];
rz(5.04866790771484) q[1];
sx q[1];
rz(3.7430735548311) q[1];
sx q[1];
rz(8.24223933219119) q[1];
rz(0.153041332960129) q[3];
sx q[3];
rz(5.85668531258638) q[3];
sx q[3];
rz(12.151910519592) q[3];
cx q[3],q[2];
rz(-2.4509584903717) q[2];
sx q[2];
rz(0.758141907053538) q[2];
sx q[2];
rz(9.49639734476014) q[2];
rz(3.2803590297699) q[3];
sx q[3];
rz(4.90571144421632) q[3];
sx q[3];
rz(11.9878566026609) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.111076809465885) q[0];
sx q[0];
rz(7.36581579049165) q[0];
sx q[0];
rz(10.9327109813611) q[0];
rz(1.12084472179413) q[1];
sx q[1];
rz(4.14506009419496) q[1];
sx q[1];
rz(6.49106166361972) q[1];
cx q[1],q[0];
rz(0.916526019573212) q[0];
sx q[0];
rz(4.34737780888612) q[0];
sx q[0];
rz(11.2367373466413) q[0];
rz(-3.29895734786987) q[2];
sx q[2];
rz(4.45562330086763) q[2];
sx q[2];
rz(10.1097751617353) q[2];
cx q[2],q[1];
rz(-1.12855911254883) q[1];
sx q[1];
rz(6.67936054070527) q[1];
sx q[1];
rz(9.96063760518237) q[1];
rz(-0.221738144755363) q[3];
sx q[3];
rz(3.37802157004411) q[3];
sx q[3];
rz(6.27234909533664) q[3];
cx q[3],q[2];
rz(3.46011090278625) q[2];
sx q[2];
rz(4.21092370350892) q[2];
sx q[2];
rz(9.92497972249194) q[2];
rz(-2.92531657218933) q[3];
sx q[3];
rz(4.62190821965272) q[3];
sx q[3];
rz(10.5342398643415) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0886539295315742) q[0];
sx q[0];
rz(1.56819883187348) q[0];
sx q[0];
rz(11.4055347204129) q[0];
rz(-0.243074640631676) q[1];
sx q[1];
rz(7.94156280358369) q[1];
sx q[1];
rz(14.5542845487516) q[1];
cx q[1],q[0];
rz(0.506211578845978) q[0];
sx q[0];
rz(2.67369091709191) q[0];
sx q[0];
rz(12.0653273820798) q[0];
rz(-4.18961715698242) q[2];
sx q[2];
rz(2.66183936794335) q[2];
sx q[2];
rz(11.843405699722) q[2];
cx q[2],q[1];
rz(-1.85503137111664) q[1];
sx q[1];
rz(1.55305448372895) q[1];
sx q[1];
rz(9.37243534847304) q[1];
rz(1.51933920383453) q[3];
sx q[3];
rz(4.11215636332566) q[3];
sx q[3];
rz(7.63320455550357) q[3];
cx q[3],q[2];
rz(2.20051765441895) q[2];
sx q[2];
rz(4.0198382457071) q[2];
sx q[2];
rz(10.2107248663823) q[2];
rz(0.032622966915369) q[3];
sx q[3];
rz(-0.918981877965383) q[3];
sx q[3];
rz(13.0552868604581) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.27633738517761) q[0];
sx q[0];
rz(5.72819724877412) q[0];
sx q[0];
rz(11.7547368764798) q[0];
rz(3.59559750556946) q[1];
sx q[1];
rz(1.78945270379121) q[1];
sx q[1];
rz(8.9586043715398) q[1];
cx q[1],q[0];
rz(2.356285572052) q[0];
sx q[0];
rz(7.94839492638642) q[0];
sx q[0];
rz(8.94027230738803) q[0];
rz(2.41455173492432) q[2];
sx q[2];
rz(5.68251458008821) q[2];
sx q[2];
rz(12.1664862394254) q[2];
cx q[2],q[1];
rz(2.58909034729004) q[1];
sx q[1];
rz(1.21181658108766) q[1];
sx q[1];
rz(7.43968710898563) q[1];
rz(3.85922265052795) q[3];
sx q[3];
rz(4.07415107091004) q[3];
sx q[3];
rz(9.98071995972797) q[3];
cx q[3],q[2];
rz(1.06352818012238) q[2];
sx q[2];
rz(1.96432760556275) q[2];
sx q[2];
rz(10.6533590316693) q[2];
rz(0.560339629650116) q[3];
sx q[3];
rz(4.4232269843393) q[3];
sx q[3];
rz(11.2413490772168) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.23357701301575) q[0];
sx q[0];
rz(2.94867393572862) q[0];
sx q[0];
rz(12.2485680341642) q[0];
rz(-0.692151486873627) q[1];
sx q[1];
rz(5.22140494187409) q[1];
sx q[1];
rz(11.6171555280606) q[1];
cx q[1],q[0];
rz(0.54587060213089) q[0];
sx q[0];
rz(3.57575795252854) q[0];
sx q[0];
rz(9.54906017928525) q[0];
rz(-2.77325415611267) q[2];
sx q[2];
rz(5.59284368355805) q[2];
sx q[2];
rz(6.07619354724094) q[2];
cx q[2],q[1];
rz(-0.918649554252625) q[1];
sx q[1];
rz(5.77841225464875) q[1];
sx q[1];
rz(10.8927836179654) q[1];
rz(3.24167823791504) q[3];
sx q[3];
rz(1.51199415524537) q[3];
sx q[3];
rz(8.19383428095981) q[3];
cx q[3],q[2];
rz(6.21583366394043) q[2];
sx q[2];
rz(4.08666005929048) q[2];
sx q[2];
rz(5.42442653178378) q[2];
rz(-0.3735211789608) q[3];
sx q[3];
rz(4.90612498124177) q[3];
sx q[3];
rz(9.50515123306915) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(6.7236123085022) q[0];
sx q[0];
rz(1.64101091225679) q[0];
sx q[0];
rz(8.47852072714969) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(2.11979722976685) q[1];
sx q[1];
rz(0.608647497492381) q[1];
sx q[1];
rz(7.98319778441592) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(1.07465445995331) q[2];
sx q[2];
rz(0.544228227930613) q[2];
sx q[2];
rz(8.30017623900577) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.597964942455292) q[3];
sx q[3];
rz(3.96313843329484) q[3];
sx q[3];
rz(15.6034722089688) q[3];
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
