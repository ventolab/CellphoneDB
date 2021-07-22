from cellphonedb.src.app.app_logger import app_logger


class LocalQueryLauncher:
    def __getattribute__(self, name):
        method = object.__getattribute__(self, name)
        if hasattr(method, '__call__'):
            app_logger.info('Launching Query {}'.format(name))

        return method

    def __init__(self, cellphonedb_app):

        self.cellphonedb_app = cellphonedb_app

    def find_interactions_by_element(self, element: str, output:str = None) -> None:
        """Search interactions in database

        This methdos recieves a list of elements as a single string separated
        by commas and uses it as a filter to search for interactions in the 
        database.

        Parameters
        ----------
        element : str
            String of elements separated by comma. 
        output : str, optional
            If provided this values is used as a path to save the query results.

        Returns
        -------
        None
            This method either outputs to the command line (if output is not 
            provided) or to a file but does not have a return value.
        
        """
        csv = self.cellphonedb_app.query.find_interactions_by_element(element).to_csv(output, index=False)
        if output:
            app_logger.info("Saving query results as '{}'".format(output))
        else:
            print(csv)

    def get_interaction_gene(self, columns: str) -> None:
        if columns:
            columns = columns.split(',')

        print(self.cellphonedb_app.query.get_all_genes(columns).to_csv(index=False))

    def autocomplete_element(self, partial_element: str) -> None:
        print(self.cellphonedb_app.query.autocomplete_launcher(partial_element).to_csv(index=False))
