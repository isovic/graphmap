// #include "make_transcript.h"
#include "index_spaced_hash_fast.h"

int IndexSpacedHashFast::MakeTranscript_(std::string annotations_path, const SequenceFile &references, SequenceFile &transcripts) {
/*
Kratke upute za koristenje SequenceFile-a.
SequenceFile drzi sve informacije o sekvencama ucitanim iz neke datoteke, moze automatski parsati FASTA, FASTQ, SAM i GFA datoteke, i to jednostavno preko konstruktora:
	SequenceFile fastqfile("sample-data/test.fastq");
	fastqfile.Verbose(stdout);	// Ispisuje neke osnovne informacije, kao npr. broj sekvenci u file-u.

Ako je potrebno proiterirati kroz sve sekvence, to se moze napraviti ovako:
	for (int64_t i=0; i<fastqfile.get_sequences().size(); i++) {
	}

Getter fastqfile.get_sequences() vraca referencu na sequences_ member koji je std::vector<SingleSequence>.
SingleSequence drzi podatke o jednoj jedinoj sekvenci (npr. u FASTA fileu koji ima vise sekvenci, ovo bi bila jedna sekvenca i jedan header; u SAM fileu, ovo bi bila jedna sekvenca u jednom retku ali takodjer ima i informacije o alibnmentu (sto u ovom trenutku nije bitno)).
Sama sekvenca (niz charactera ACTG baza) spremljena je u obliku int8_t* data_ member polja, pristupas mu preko .get_data() gettera. Ovo polje mozes cast-ati direktno u (char *) ako je potrebno, ali nije \0 terminiran. Velicinu data polja dobijes s .get_data_length()

	for (int64_t i=0; i<fastqfile.get_sequences().size(); i++) {
		auto seq = fastqfile.get_sequences()[i];	// Pointer na SingleSequence objekt.

		// Ispisivanje pojedinih baza u sekvenci.
		for (int64_t j=0; j<seq->get_data_length(); j++) {
			printf ("%c", (char) seq->get_data()[j]);
		}
		printf ("\n");
	}

Na ovaj nacin mozes pristupati svim postojecim ucitanim referentnim sekvencama, one ce biti predane u references objektu.

Kada generiras jedan transkriptom, dodaj ga u transcripts objekt na sljedeci nacin.
Recimo da imas neki string koji sadrzi tvoju sekvencu (cisto radi jednostavnosti primjera):
	std::string sekvenca("ACTGCTGCTA");
	std::string header("Sekvenca 1");

Prvo, generiraj novi SingleSequence objekt:
	SingleSequence *s = new SingleSequence;
	s->InitHeaderAndDataFromAscii(header.c_str(), header.size(), sekvenca.c_str(), sekvenca.size(), 1);	// Zadnji broj 1 je proizvoljni ID sekvence, ali bilo bi dobro da je po redu.

Dodati taj objekt u transcripts:
	transcripts.AddSequence(s, true);	// "true" kaze SequenceFile objektu da kod destrukcije oslobodi memoriju od "s" objekta.


/////////////////////////////////////
/////////////////////////////////////

Plan je ovakav: IndexSpacedHashFast ce ostati skoro isti kao do sada s jednom razlikom, ako je commanline parametar specificirao da ce na ulazu biti i anotirani file s exonima (GFF ili BED), pozvat cu prvo ovu funkciju.
SequenceFile transcripts, koji ces popuniti ovdje, cu ja onda indeksirati umjesto pravih referenci.
Slobodno dodaj membere u Index koji ti trebaju, jer ces morati imati dodatne info kako bi mogao napraviti konverziju nazad iz transcriptome prostora u prostor genoma.
Takodjer, slobodno dodaj parametre u sucelju gore, s obzirom da dok ne implementiras inicijalnu verziju ne mogu niti ja implementirati u kod, sve je jos u zraku :-)

S obzirom da moras implementirati i konverziju iz koordinata alignmenta nazad u transkriptom, taj dio je trenutno jos malo tricky. U GraphMap-u imam strukturu podataka koja drzi cijeli alignment, ali moram opisati tocno kako da je koristis.
Dok ja to ne napravim, ti u medjuvremenu implementiraj ovu funkcionalnost.

Ako bilo sto od ovoga iznad nije jasno, samo pitaj.

S obzirom da vec imas rijesenu implementaciju u Pythonu, mislim da ovdje ne bi trebalo biti previse posla, bilo bi odlicno kada bi to uspjeli finalizirati cim prije i staviti ljudima na koristenje.

*/

	transcripts.Clear();



	return 0;
}
