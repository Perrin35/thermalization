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
rz(0.736833572387695) q[0];
sx q[0];
rz(4.92173043091828) q[0];
sx q[0];
rz(11.1877274274747) q[0];
rz(2.28400754928589) q[1];
sx q[1];
rz(4.62559142907197) q[1];
sx q[1];
rz(8.97388226389095) q[1];
cx q[1],q[0];
rz(4.88440561294556) q[0];
sx q[0];
rz(3.23158116837079) q[0];
sx q[0];
rz(11.159801220886) q[0];
rz(-1.45920431613922) q[2];
sx q[2];
rz(1.09422317345674) q[2];
sx q[2];
rz(6.28056452273532) q[2];
cx q[2],q[1];
rz(3.29708695411682) q[1];
sx q[1];
rz(4.93440738518769) q[1];
sx q[1];
rz(7.34199211596652) q[1];
rz(5.59941816329956) q[3];
sx q[3];
rz(4.78658083279664) q[3];
sx q[3];
rz(8.56578061579868) q[3];
cx q[3],q[2];
rz(-1.98400545120239) q[2];
sx q[2];
rz(4.82365134556825) q[2];
sx q[2];
rz(10.2690666079442) q[2];
rz(-0.441011607646942) q[3];
sx q[3];
rz(3.49725812871987) q[3];
sx q[3];
rz(8.81875458954974) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.592502951622009) q[0];
sx q[0];
rz(4.37140849431092) q[0];
sx q[0];
rz(9.68786979316875) q[0];
rz(-0.943537652492523) q[1];
sx q[1];
rz(2.54480055172975) q[1];
sx q[1];
rz(10.611010169975) q[1];
cx q[1],q[0];
rz(-1.10375583171844) q[0];
sx q[0];
rz(3.08259664301807) q[0];
sx q[0];
rz(8.17763302325412) q[0];
rz(-0.0621210597455502) q[2];
sx q[2];
rz(4.51336184342439) q[2];
sx q[2];
rz(11.2635206937711) q[2];
cx q[2],q[1];
rz(-3.40304493904114) q[1];
sx q[1];
rz(4.83746174176271) q[1];
sx q[1];
rz(11.8997788190763) q[1];
rz(-0.192477867007256) q[3];
sx q[3];
rz(5.39821091492707) q[3];
sx q[3];
rz(4.82610032557651) q[3];
cx q[3],q[2];
rz(-1.01203227043152) q[2];
sx q[2];
rz(2.13887813885743) q[2];
sx q[2];
rz(13.7258400678556) q[2];
rz(-3.51267743110657) q[3];
sx q[3];
rz(4.77871218522126) q[3];
sx q[3];
rz(9.73571387528583) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.761125802993774) q[0];
sx q[0];
rz(2.01158705552156) q[0];
sx q[0];
rz(7.08990643023654) q[0];
rz(-0.21356788277626) q[1];
sx q[1];
rz(2.64533058007295) q[1];
sx q[1];
rz(8.60455832480594) q[1];
cx q[1],q[0];
rz(-1.85042822360992) q[0];
sx q[0];
rz(3.83680007060105) q[0];
sx q[0];
rz(9.24887797831699) q[0];
rz(4.49325752258301) q[2];
sx q[2];
rz(4.22886982758576) q[2];
sx q[2];
rz(8.90371272563144) q[2];
cx q[2],q[1];
rz(0.0809382349252701) q[1];
sx q[1];
rz(5.06456175644929) q[1];
sx q[1];
rz(8.97459474801227) q[1];
rz(1.68882882595062) q[3];
sx q[3];
rz(7.391994150477) q[3];
sx q[3];
rz(9.42997893354996) q[3];
cx q[3],q[2];
rz(-0.310725808143616) q[2];
sx q[2];
rz(4.64222160180146) q[2];
sx q[2];
rz(13.4971704244535) q[2];
rz(-0.155497744679451) q[3];
sx q[3];
rz(4.77953139145906) q[3];
sx q[3];
rz(6.57473728655978) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.75472331047058) q[0];
sx q[0];
rz(3.86297211249406) q[0];
sx q[0];
rz(10.3360531091611) q[0];
rz(-3.57991194725037) q[1];
sx q[1];
rz(4.46378353436524) q[1];
sx q[1];
rz(10.7452029943387) q[1];
cx q[1],q[0];
rz(-2.51169180870056) q[0];
sx q[0];
rz(1.98038473923738) q[0];
sx q[0];
rz(10.9843125104825) q[0];
rz(0.206409737467766) q[2];
sx q[2];
rz(5.1552461703592) q[2];
sx q[2];
rz(10.129539525501) q[2];
cx q[2],q[1];
rz(-0.988523244857788) q[1];
sx q[1];
rz(1.30779019196565) q[1];
sx q[1];
rz(10.2911045312802) q[1];
rz(-2.35841989517212) q[3];
sx q[3];
rz(5.16765299637849) q[3];
sx q[3];
rz(10.0406999349515) q[3];
cx q[3],q[2];
rz(-2.03582549095154) q[2];
sx q[2];
rz(2.21290090878541) q[2];
sx q[2];
rz(9.08239684104129) q[2];
rz(0.176774471998215) q[3];
sx q[3];
rz(6.71632209618623) q[3];
sx q[3];
rz(7.42415831088229) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.311556100845337) q[0];
sx q[0];
rz(0.727685840921946) q[0];
sx q[0];
rz(11.7010800599973) q[0];
rz(4.36814165115356) q[1];
sx q[1];
rz(5.29391613801057) q[1];
sx q[1];
rz(13.8669872045438) q[1];
cx q[1],q[0];
rz(-1.96017014980316) q[0];
sx q[0];
rz(1.92480650742585) q[0];
sx q[0];
rz(9.39788896813198) q[0];
rz(-1.9373539686203) q[2];
sx q[2];
rz(1.50368109543855) q[2];
sx q[2];
rz(10.8221217155378) q[2];
cx q[2],q[1];
rz(1.27582669258118) q[1];
sx q[1];
rz(7.52566543419892) q[1];
sx q[1];
rz(7.28396055697604) q[1];
rz(-2.99556136131287) q[3];
sx q[3];
rz(4.9720049222284) q[3];
sx q[3];
rz(9.12654463051959) q[3];
cx q[3],q[2];
rz(-2.13179230690002) q[2];
sx q[2];
rz(4.29100719292695) q[2];
sx q[2];
rz(13.0435702562253) q[2];
rz(-0.192084327340126) q[3];
sx q[3];
rz(4.83527830441529) q[3];
sx q[3];
rz(8.49166128634616) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.796482980251312) q[0];
sx q[0];
rz(6.89744821389253) q[0];
sx q[0];
rz(9.43652826211556) q[0];
rz(0.55039644241333) q[1];
sx q[1];
rz(4.49791375001008) q[1];
sx q[1];
rz(7.83633503913089) q[1];
cx q[1],q[0];
rz(-3.30045485496521) q[0];
sx q[0];
rz(5.04749885399873) q[0];
sx q[0];
rz(5.91085646151706) q[0];
rz(-1.59266519546509) q[2];
sx q[2];
rz(4.34658232529695) q[2];
sx q[2];
rz(7.36260817050144) q[2];
cx q[2],q[1];
rz(2.46588182449341) q[1];
sx q[1];
rz(4.79721084435517) q[1];
sx q[1];
rz(5.42243859767123) q[1];
rz(-0.0879004150629044) q[3];
sx q[3];
rz(3.84631166060502) q[3];
sx q[3];
rz(12.2288694143216) q[3];
cx q[3],q[2];
rz(-2.63402962684631) q[2];
sx q[2];
rz(2.48122486670549) q[2];
sx q[2];
rz(8.14222679137393) q[2];
rz(1.36986184120178) q[3];
sx q[3];
rz(4.88785007794435) q[3];
sx q[3];
rz(10.5432348012845) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.581055283546448) q[0];
sx q[0];
rz(3.30955437024171) q[0];
sx q[0];
rz(10.1020375251691) q[0];
rz(3.29340028762817) q[1];
sx q[1];
rz(1.37445238431031) q[1];
sx q[1];
rz(11.5893213510434) q[1];
cx q[1],q[0];
rz(2.55172944068909) q[0];
sx q[0];
rz(4.78771200974519) q[0];
sx q[0];
rz(7.26079294680759) q[0];
rz(1.57004654407501) q[2];
sx q[2];
rz(3.27363348205621) q[2];
sx q[2];
rz(9.537206670634) q[2];
cx q[2],q[1];
rz(-0.207714200019836) q[1];
sx q[1];
rz(1.90954509575898) q[1];
sx q[1];
rz(7.77229306697055) q[1];
rz(0.81959593296051) q[3];
sx q[3];
rz(3.59573808510835) q[3];
sx q[3];
rz(8.79608646630451) q[3];
cx q[3],q[2];
rz(-1.38924217224121) q[2];
sx q[2];
rz(3.96328738530213) q[2];
sx q[2];
rz(8.41205224990054) q[2];
rz(-1.18794536590576) q[3];
sx q[3];
rz(4.21411159833009) q[3];
sx q[3];
rz(8.93756592868968) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.12412071228027) q[0];
sx q[0];
rz(6.24982467492158) q[0];
sx q[0];
rz(10.1234377980153) q[0];
rz(1.12208044528961) q[1];
sx q[1];
rz(3.98768499692018) q[1];
sx q[1];
rz(8.17539236544772) q[1];
cx q[1],q[0];
rz(-0.699853360652924) q[0];
sx q[0];
rz(6.06732860405976) q[0];
sx q[0];
rz(8.07178995608493) q[0];
rz(0.72624260187149) q[2];
sx q[2];
rz(5.62659064133699) q[2];
sx q[2];
rz(9.51114259510442) q[2];
cx q[2],q[1];
rz(0.577512085437775) q[1];
sx q[1];
rz(3.11875935097272) q[1];
sx q[1];
rz(8.09900579451724) q[1];
rz(0.384371548891068) q[3];
sx q[3];
rz(0.669251831369944) q[3];
sx q[3];
rz(8.59007606505557) q[3];
cx q[3],q[2];
rz(-0.329687833786011) q[2];
sx q[2];
rz(5.49702087243135) q[2];
sx q[2];
rz(10.6032720565717) q[2];
rz(1.45684361457825) q[3];
sx q[3];
rz(5.22075024445588) q[3];
sx q[3];
rz(9.8069131731908) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.64173984527588) q[0];
sx q[0];
rz(4.52116707165773) q[0];
sx q[0];
rz(10.7177903413694) q[0];
rz(-1.71998429298401) q[1];
sx q[1];
rz(4.17797449429566) q[1];
sx q[1];
rz(11.9687928914945) q[1];
cx q[1],q[0];
rz(2.98747277259827) q[0];
sx q[0];
rz(6.14264622529084) q[0];
sx q[0];
rz(6.56181809901401) q[0];
rz(-0.525528371334076) q[2];
sx q[2];
rz(4.21843674977357) q[2];
sx q[2];
rz(10.675688481323) q[2];
cx q[2],q[1];
rz(-0.240391954779625) q[1];
sx q[1];
rz(1.24820926983888) q[1];
sx q[1];
rz(7.54028627871677) q[1];
rz(1.68652760982513) q[3];
sx q[3];
rz(5.69413343270356) q[3];
sx q[3];
rz(11.655866599075) q[3];
cx q[3],q[2];
rz(0.222758159041405) q[2];
sx q[2];
rz(4.83178046544129) q[2];
sx q[2];
rz(8.1914180278699) q[2];
rz(5.38179922103882) q[3];
sx q[3];
rz(3.26165025879676) q[3];
sx q[3];
rz(7.78146216868564) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.0478865876793861) q[0];
sx q[0];
rz(3.91355028946931) q[0];
sx q[0];
rz(9.44843821264013) q[0];
rz(0.95611447095871) q[1];
sx q[1];
rz(4.45119098027284) q[1];
sx q[1];
rz(11.894198870651) q[1];
cx q[1],q[0];
rz(-1.72608721256256) q[0];
sx q[0];
rz(3.1531102844798) q[0];
sx q[0];
rz(9.90175337194606) q[0];
rz(0.00928795430809259) q[2];
sx q[2];
rz(2.08229676087434) q[2];
sx q[2];
rz(8.01952753066226) q[2];
cx q[2],q[1];
rz(-0.0862395539879799) q[1];
sx q[1];
rz(1.90859106381471) q[1];
sx q[1];
rz(10.2783196926038) q[1];
rz(-2.00362944602966) q[3];
sx q[3];
rz(1.43086174328858) q[3];
sx q[3];
rz(9.82111541032001) q[3];
cx q[3],q[2];
rz(-0.512229263782501) q[2];
sx q[2];
rz(5.05625608761842) q[2];
sx q[2];
rz(9.79473657011195) q[2];
rz(1.63791787624359) q[3];
sx q[3];
rz(4.02748558123643) q[3];
sx q[3];
rz(11.3653987407605) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.17553627490997) q[0];
sx q[0];
rz(4.37310984929139) q[0];
sx q[0];
rz(11.1302566289823) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(2.90803408622742) q[1];
sx q[1];
rz(5.42954746087129) q[1];
sx q[1];
rz(7.1030630826871) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-0.313793867826462) q[2];
sx q[2];
rz(7.50068727334077) q[2];
sx q[2];
rz(9.00291634201213) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(3.2831666469574) q[3];
sx q[3];
rz(5.15133014519746) q[3];
sx q[3];
rz(10.7782805919568) q[3];
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